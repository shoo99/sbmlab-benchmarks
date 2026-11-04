## Builds the doublet-detection inputs with their truth from the downloads (issue #7).
##
## Each data set goes to <work_dir>/<name>/:
##   counts.mtx.gz, barcodes.tsv, genes.tsv   genes x droplets, every droplet of the filtered matrix
##   truth.tsv                                barcode, truth (doublet / singlet / unknown)
##   info.tsv                                 droplets, expected doublet rate, truth counts, what the truth can see
## Data sets and truth
##   kang_ctrl, kang_stim   Kang et al. 2018 (GSE96583 batch 2), 8 donors pooled per lane; demuxlet calls droplets with
##                          two genotypes "doublet" (truth doublet), "singlet" (truth singlet), "ambs" (unknown).
##                          Doublets of two cells from one donor are invisible to this truth: share 1 - sum(p_donor^2)
##                          of all doublets is visible.
##   hgmm_12k               10x Genomics human (HEK293T) + mouse (NIH3T3) mixture, Cell Ranger 2.1 filtered matrices.
##                          Minor-species share of UMIs >= 0.2: cross-species doublet; <= 0.1: singlet; else unknown
##                          (ambient RNA gives a few percent). Same-species doublets are invisible: about half of all.
## Expected doublet rate: 0.8% per 1,000 droplets (the 10x rule of thumb), capped at 25%; EXPECTED_RATE overrides it.
##
## Usage: Rscript prepare.R <data_dir> <work_dir> [dataset ...]   (default: kang_ctrl kang_stim hgmm_12k)
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: Rscript prepare.R <data_dir> <work_dir> [dataset ...]")
data_dir <- args[1]; work_dir <- args[2]
datasets <- if (length(args) > 2) args[-(1:2)] else c("kang_ctrl", "kang_stim", "hgmm_12k")
MINOR_DOUBLET <- 0.2; MINOR_SINGLET <- 0.1

expected_rate <- function(n) if (nzchar(Sys.getenv("EXPECTED_RATE"))) as.numeric(Sys.getenv("EXPECTED_RATE")) else min(0.25, 0.008 * n / 1000)
read_lines_any <- function(path) if (grepl("\\.gz$", path)) readLines(gzfile(path)) else readLines(path)
first_existing <- function(dir, names) { p <- file.path(dir, names); p <- p[file.exists(p)]; if (length(p)) p[1] else NA }

write_dataset <- function(name, counts, truth, visible, note) {
  d <- file.path(work_dir, name); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(d, "counts.mtx")
  invisible(writeMM(counts, mtx))
  if (file.exists(paste0(mtx, ".gz"))) file.remove(paste0(mtx, ".gz"))
  system2("gzip", c("-n", shQuote(mtx)))
  writeLines(colnames(counts), file.path(d, "barcodes.tsv"))
  writeLines(rownames(counts), file.path(d, "genes.tsv"))
  utils::write.table(data.frame(barcode = colnames(counts), truth = truth), file.path(d, "truth.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  info <- data.frame(dataset = name, droplets = ncol(counts), genes = nrow(counts),
                     expected_rate = round(expected_rate(ncol(counts)), 4), truth_doublet = sum(truth == "doublet"),
                     truth_singlet = sum(truth == "singlet"), truth_unknown = sum(truth == "unknown"),
                     visible_share = round(visible, 3), truth_note = note)
  utils::write.table(info, file.path(d, "info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-10s %6d droplets, %5d genes; truth: %d doublets, %d singlets, %d unknown; expected rate %.3f; truth sees %.0f%% of doublets\n",
              name, ncol(counts), nrow(counts), info$truth_doublet, info$truth_singlet, info$truth_unknown,
              info$expected_rate, 100 * visible))
}

## ---- Kang 2018: demuxlet labels
prepare_kang <- function(cond) {
  files <- list(ctrl = c("GSM2560248_2.1.mtx.gz", "GSM2560248_barcodes.tsv.gz"),
                stim = c("GSM2560249_2.2.mtx.gz", "GSM2560249_barcodes.tsv.gz"))[[cond]]
  counts <- as(readMM(gzfile(file.path(data_dir, files[1]))), "CsparseMatrix")
  barcodes <- readLines(gzfile(file.path(data_dir, files[2])))
  genes <- utils::read.delim(gzfile(file.path(data_dir, "GSE96583_batch2.genes.tsv.gz")), header = FALSE)
  stopifnot(nrow(counts) == nrow(genes), ncol(counts) == length(barcodes))
  dimnames(counts) <- list(make.unique(genes[[1]]), barcodes)
  cd <- utils::read.delim(gzfile(file.path(data_dir, "GSE96583_batch2.total.tsne.df.tsv.gz")), stringsAsFactors = FALSE)
  if (names(cd)[1] %in% c("X", "barcode", "cell_barcode")) { rownames(cd) <- cd[[1]]; cd <- cd[, -1] }
  cd <- cd[cd$stim == cond, ]
  idx <- match(barcodes, rownames(cd))
  if (mean(!is.na(idx)) < 0.8) idx <- match(sub("-[0-9]+$", "", barcodes), sub("-[0-9]+$", "", rownames(cd)))
  if (mean(!is.na(idx)) < 0.8) stop(sprintf("kang_%s: only %.0f%% of barcodes match the demuxlet table", cond, 100 * mean(!is.na(idx))))
  call <- cd$multiplets[idx]
  truth <- ifelse(is.na(call), "unknown", ifelse(call == "doublet", "doublet", ifelse(call == "singlet", "singlet", "unknown")))
  p <- prop.table(table(cd$ind[idx][truth == "singlet"]))
  write_dataset(paste0("kang_", cond), counts, truth, 1 - sum(p^2),
                "demuxlet: droplets with two genotypes; same-donor doublets invisible")
}

## ---- 10x human + mouse mixture: species share of UMIs
read_10x <- function(dir) {
  mtx <- first_existing(dir, c("matrix.mtx.gz", "matrix.mtx"))
  feat <- first_existing(dir, c("features.tsv.gz", "genes.tsv.gz", "genes.tsv", "features.tsv"))
  bc <- first_existing(dir, c("barcodes.tsv.gz", "barcodes.tsv"))
  m <- as(readMM(if (grepl("gz$", mtx)) gzfile(mtx) else mtx), "CsparseMatrix")
  g <- utils::read.delim(if (grepl("gz$", feat)) gzfile(feat) else feat, header = FALSE)
  dimnames(m) <- list(g[[1]], read_lines_any(bc))
  m
}
prepare_hgmm <- function(name) {
  base <- file.path(data_dir, name, "filtered_gene_bc_matrices")
  if (dir.exists(base)) {                    # Cell Ranger 2: one folder per genome
    parts <- lapply(c(human = "hg19", mouse = "mm10"), function(g) read_10x(file.path(base, g)))
    bcs <- sort(unique(unlist(lapply(parts, colnames))))
    parts <- lapply(names(parts), function(sp) {
      m <- parts[[sp]]
      full <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(nrow(m), length(bcs)),
                           dimnames = list(paste0(sp, "_", rownames(m)), bcs))
      full[, colnames(m)] <- m
      full
    })
    counts <- rbind(parts[[1]], parts[[2]])
  } else {                                   # Cell Ranger 3+: one matrix, gene ids prefixed with the genome
    counts <- read_10x(file.path(data_dir, name, "filtered_feature_bc_matrix"))
    rownames(counts) <- sub("^(hg19|GRCh38)_+", "human_", sub("^(mm10|GRCm38|GRCm39)_+", "mouse_", rownames(counts)))
  }
  human <- Matrix::colSums(counts[startsWith(rownames(counts), "human_"), , drop = FALSE])
  mouse <- Matrix::colSums(counts[startsWith(rownames(counts), "mouse_"), , drop = FALSE])
  minor <- pmin(human, mouse) / pmax(human + mouse, 1)
  truth <- ifelse(minor >= MINOR_DOUBLET, "doublet", ifelse(minor <= MINOR_SINGLET, "singlet", "unknown"))
  ph <- mean(human[truth == "singlet"] > mouse[truth == "singlet"])
  write_dataset(name, counts, truth, 2 * ph * (1 - ph),
                sprintf("species share of UMIs: minor >= %.2f doublet, <= %.2f singlet; same-species doublets invisible",
                        MINOR_DOUBLET, MINOR_SINGLET))
}

for (ds in datasets) {
  if (ds %in% c("kang_ctrl", "kang_stim")) prepare_kang(sub("kang_", "", ds)) else prepare_hgmm(ds)
}
cat(R.version.string, " | Matrix ", format(packageVersion("Matrix")), "\n", sep = "")
