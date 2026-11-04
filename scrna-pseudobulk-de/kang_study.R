## Cell-level vs pseudobulk DE on real PBMC from 8 donors: Kang et al. 2018 control cells (issue #6, part 2)
##
## Kang et al. (2018, Nat Biotechnol 36:89; GEO GSE96583, batch 2) pooled PBMC from 8 lupus patients and
## demultiplexed them by genotype (demuxlet), so every cell has a donor and a doublet call.
## Only the unstimulated (control) cells are used, so the truth is set by the design:
##   null     the 8 donors are split 4 vs 4; all 35 distinct splits are run. Nothing separates the groups except
##            which people are in them, so every call is false.
##   planted  for N_PLANTED of the splits, 10% of the tested genes are changed 2-fold by binomial thinning: counts
##            in group A cells are thinned to 1/2 for half of them (B higher), counts in group B cells for the
##            other half (B lower). Thinning keeps the count nature of the data and changes only the chosen genes.
## Cells: demuxlet singlets with a cell-type label; cell types with at least MIN_CELLS cells in every donor.
## Genes: detected in at least 5% of the cell type's control cells; all methods test the same genes.
##
## Usage: Rscript kang_study.R [data_dir] > kang_study.out    (data_dir: what download.sh fetched; default data/)
##   env: N_SPLITS (35), N_PLANTED (10), MIN_CELLS (20), KANG_CORES (1), DE_MAST=1 to add MAST, OUT_DIR (this folder)
## Writes kang_metrics.csv (one row per cell type, split, scenario and method); score.R adds it to scores.csv.
## Seed: set.seed(2018) before the planted genes of each cell type are drawn; thinning uses set.seed(2018 + split).
suppressPackageStartupMessages({ library(Matrix); library(parallel) })
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "de_methods.R"))
options(width = 140)
args <- commandArgs(TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path(script_dir, "data")
N_SPLITS <- as.integer(Sys.getenv("N_SPLITS", "35"))
N_PLANTED <- as.integer(Sys.getenv("N_PLANTED", "10"))
MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "20"))
CORES <- as.integer(Sys.getenv("KANG_CORES", "1"))
OUT_DIR <- Sys.getenv("OUT_DIR", script_dir)
FC <- 2

## ---- read the control matrix and the demuxlet / cell-type table
f <- function(name) file.path(data_dir, name)
counts <- as(readMM(gzfile(f("GSM2560248_2.1.mtx.gz"))), "CsparseMatrix")
genes <- utils::read.delim(gzfile(f("GSE96583_batch2.genes.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
barcodes <- readLines(gzfile(f("GSM2560248_barcodes.tsv.gz")))
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == length(barcodes))
rownames(counts) <- genes[[1]]; colnames(counts) <- barcodes
cd <- utils::read.delim(gzfile(f("GSE96583_batch2.total.tsne.df.tsv.gz")), stringsAsFactors = FALSE)
if (!"ind" %in% names(cd)) stop("the cell table has no 'ind' (donor) column")
if (names(cd)[1] %in% c("X", "barcode", "cell_barcode")) { rownames(cd) <- cd[[1]]; cd <- cd[, -1] }  # barcodes as a column
cd <- cd[cd$stim == "ctrl", ]
## the table's row names are barcodes; match them exactly, or after dropping a trailing "-<digits>" on both sides
idx <- match(barcodes, rownames(cd))
if (mean(!is.na(idx)) < 0.8) idx <- match(sub("-[0-9]+$", "", barcodes), sub("-[0-9]+$", "", rownames(cd)))
cat(sprintf("== Kang 2018 control: %d cells in the matrix, %d matched to the demuxlet / cell-type table (%.1f%%)\n",
            length(barcodes), sum(!is.na(idx)), 100 * mean(!is.na(idx))))
if (mean(!is.na(idx)) < 0.8) stop("fewer than 80% of the barcodes match the table; check the downloaded files")
meta <- cd[idx, c("ind", "cell", "multiplets")]
cat("   demuxlet calls:", paste(names(table(meta$multiplets, useNA = "ifany")), table(meta$multiplets, useNA = "ifany"),
                               sep = " ", collapse = ", "), "\n")
ok <- !is.na(meta$multiplets) & meta$multiplets == "singlet" & !is.na(meta$cell) & meta$cell != ""
counts <- counts[, ok]; meta <- meta[ok, ]
meta$donor <- factor(meta$ind)
donors <- levels(meta$donor)
stopifnot(length(donors) == 8)
per <- table(meta$cell, meta$donor)
types <- rownames(per)[apply(per, 1, min) >= MIN_CELLS]
cat(sprintf("   singlets: %d in %d donors; cell types with >= %d cells in every donor: %s\n", ncol(counts), length(donors),
            MIN_CELLS, paste(sprintf("%s (%d)", types, rowSums(per[types, , drop = FALSE])), collapse = ", ")))

## ---- all 35 balanced splits of 8 donors (the first donor always in group A)
splits <- combn(donors[-1], 3, function(x) c(donors[1], x), simplify = FALSE)
splits <- splits[seq_len(min(N_SPLITS, length(splits)))]
set.seed(2018)
planted_splits <- sort(sample(seq_along(splits), min(N_PLANTED, length(splits))))

thin <- function(x, rows, cols, p) {          # binomial thinning of x[rows, cols] to share p
  sub <- x[rows, cols, drop = FALSE]
  sub@x <- as.numeric(stats::rbinom(length(sub@x), sub@x, p))
  x[rows, cols] <- sub
  drop0(x)
}

out <- list()
for (ct in types) {
  cells <- meta$cell == ct
  x <- counts[, cells]; donor <- droplevels(meta$donor[cells])
  x <- x[Matrix::rowMeans(x > 0) >= 0.05, ]
  set.seed(2018)
  pl <- sample(nrow(x), round(nrow(x) / 10))
  up <- pl[seq_len(length(pl) %/% 2)]; down <- setdiff(pl, up)
  truth <- seq_len(nrow(x)) %in% pl
  t0 <- proc.time()[[3]]
  jobs <- c(lapply(seq_along(splits), function(s) list(s = s, planted = FALSE)),
            lapply(planted_splits, function(s) list(s = s, planted = TRUE)))
  res <- mclapply(jobs, function(j) {
    group <- factor(ifelse(donor %in% splits[[j$s]], "A", "B"), levels = c("A", "B"))
    y <- x; tr <- rep(FALSE, nrow(x))
    if (j$planted) {
      set.seed(2018 + j$s)
      y <- thin(y, up, which(group == "A"), 1 / FC)
      y <- thin(y, down, which(group == "B"), 1 / FC)
      tr <- truth
    }
    p <- run_methods(y, donor, group)
    st <- t(apply(p, 2, rep_stats, truth = tr))
    data.frame(cell_type = ct, cells = ncol(x), genes = nrow(x), scenario = if (j$planted) "planted" else "null",
               split = paste(splits[[j$s]], collapse = "+"), method = rownames(st), st, row.names = NULL)
  }, mc.cores = CORES)
  res <- do.call(rbind, res)
  out[[ct]] <- res
  s <- aggregate(cbind(p05, false, power, fdp) ~ scenario + method, data = res, FUN = mean, na.action = na.pass)
  s <- s[order(s$scenario, match(s$method, de_methods())), ]
  cat(sprintf("\n== %s: %d cells, %d genes tested; %d null splits, %d planted (%.0f s)\n", ct, ncol(x), nrow(x),
              length(splits), length(planted_splits), proc.time()[[3]] - t0))
  s[, -(1:2)] <- round(s[, -(1:2)], 3)
  print(s, row.names = FALSE)
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(OUT_DIR, "kang_metrics.csv"), row.names = FALSE)
version <- sprintf("edgeR %s; DESeq2 %s; limma %s; R %s (Wilcoxon rank-sum, Welch t)%s", packageVersion("edgeR"),
                   packageVersion("DESeq2"), packageVersion("limma"), getRversion(),
                   if (use_mast()) paste0("; MAST ", packageVersion("MAST")) else "")
utils::write.csv(data.frame(version = version), file.path(OUT_DIR, "kang_versions.csv"), row.names = FALSE)
cat(sprintf("\nwrote kang_metrics.csv: %d rows, and kang_versions.csv\n", nrow(res)))
cat(R.version.string, " | edgeR ", format(packageVersion("edgeR")), " | DESeq2 ", format(packageVersion("DESeq2")),
    " | limma ", format(packageVersion("limma")), if (use_mast()) paste(" | MAST", packageVersion("MAST")), "\n", sep = "")
