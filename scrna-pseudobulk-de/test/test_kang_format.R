## Known-answer test of kang_study.R on synthetic data in the format of the GSE96583 batch-2 files (run by CI).
##
## 8 donors x 2 cell types, 60-120 control cells per donor and type, 600 genes. Each gene has a donor effect
## (log-normal, SD 0.3), so cells of one donor are not independent replicates. The cell table also holds doublets
## and stimulated cells, which kang_study.R must leave out. Known answers:
##   null splits      cell-level Wilcoxon calls many genes (share p < 0.05 >= 0.15); edgeR and limma-voom stay
##                    near 0.05 (<= 0.10)
##   planted 2-fold   edgeR finds >= 30% of the changed genes with FDP <= 0.2
##   bookkeeping      only control singlets are used; one row per cell type, split, scenario and method
## Usage: Rscript test/test_kang_format.R   (about 1 min; writes only to a temporary directory)
suppressPackageStartupMessages(library(Matrix))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
folder <- normalizePath(file.path(script_dir, ".."))
tmp <- tempfile("synthetic-kang-"); dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
set.seed(5)
donors <- c("101", "1015", "1016", "1039", "107", "1244", "1256", "1488")
types <- c("CD14+ Monocytes", "B cells")
G <- 600
mu <- rlnorm(G, -0.5, 1.2)
type_eff <- matrix(rlnorm(G * 2, 0, 0.5), G)
cells <- do.call(rbind, lapply(donors, function(d) do.call(rbind, lapply(1:2, function(t)
  data.frame(ind = d, cell = types[t], n = sample(60:120, 1))))))
cells <- cells[rep(seq_len(nrow(cells)), cells$n), c("ind", "cell")]
n <- nrow(cells)
donor_eff <- matrix(rlnorm(G * 8, 0, 0.3), G, dimnames = list(NULL, donors))
m <- mu * type_eff[, match(cells$cell, types)] * donor_eff[, cells$ind] * rep(rlnorm(n, 0, 0.3), each = G)
counts <- Matrix(matrix(rnbinom(length(m), mu = m, size = 5), G), sparse = TRUE)
barcode <- sprintf("%s-1", replicate(n, paste(sample(c("A", "C", "G", "T"), 14, TRUE), collapse = "")))
while (anyDuplicated(barcode)) barcode[duplicated(barcode)] <- sprintf("%s-1", paste(sample(c("A", "C", "G", "T"), 14, TRUE), collapse = ""))
cells$multiplets <- ifelse(runif(n) < 0.1, "doublet", "singlet")
invisible(writeMM(counts, file.path(tmp, "GSM2560248_2.1.mtx")))
system2("gzip", file.path(tmp, "GSM2560248_2.1.mtx"))
writeLines(barcode, gzfile(file.path(tmp, "GSM2560248_barcodes.tsv.gz")))
write.table(data.frame(sprintf("ENSG%011d", seq_len(G)), sprintf("GENE%d", seq_len(G))),
            gzfile(file.path(tmp, "GSE96583_batch2.genes.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
tab <- data.frame(tsne1 = 0, tsne2 = 0, ind = cells$ind, stim = "ctrl", cluster = 1, cell = cells$cell,
                  multiplets = cells$multiplets, row.names = barcode)
stim <- data.frame(tsne1 = rep(0, 50), tsne2 = 0, ind = donors[1], stim = "stim", cluster = 1, cell = types[1],
                   multiplets = "singlet", row.names = sprintf("STIM%05d-1", 1:50))
write.table(rbind(tab, stim), gzfile(file.path(tmp, "GSE96583_batch2.total.tsne.df.tsv.gz")), sep = "\t", quote = FALSE)

log <- system2("Rscript", c(shQuote(file.path(folder, "kang_study.R")), shQuote(tmp)), stdout = TRUE, stderr = TRUE,
               env = c("N_SPLITS=6", "N_PLANTED=3", paste0("OUT_DIR=", shQuote(tmp)), "OPENBLAS_NUM_THREADS=1"))
cat(log, sep = "\n")
if (!is.null(attr(log, "status"))) stop("kang_study.R failed")
k <- read.csv(file.path(tmp, "kang_metrics.csv"), stringsAsFactors = FALSE)
v <- function(ct, m, sc, col) mean(k[k$cell_type == ct & k$method == m & k$scenario == sc, col])
singlets <- table(cells$cell[cells$multiplets == "singlet"])
checks <- list()
for (ct in types) {
  checks <- c(checks, list(
    list(sprintf("%s: null, Wilcoxon (cells) share p < 0.05 >= 0.15", ct), v(ct, "Wilcoxon (cells)", "null", "p05") >= 0.15),
    list(sprintf("%s: null, edgeR QL share p < 0.05 <= 0.10", ct), v(ct, "edgeR QL (pseudobulk)", "null", "p05") <= 0.10),
    list(sprintf("%s: null, limma-voom share p < 0.05 <= 0.10", ct), v(ct, "limma-voom (pseudobulk)", "null", "p05") <= 0.10),
    list(sprintf("%s: planted, edgeR QL power >= 0.3", ct), v(ct, "edgeR QL (pseudobulk)", "planted", "power") >= 0.3),
    list(sprintf("%s: planted, edgeR QL FDP <= 0.2", ct), v(ct, "edgeR QL (pseudobulk)", "planted", "fdp") <= 0.2),
    list(sprintf("%s: only control singlets used (%d cells)", ct, singlets[[ct]]), all(k$cells[k$cell_type == ct] == singlets[[ct]]))))
}
checks <- c(checks, list(list("one row per cell type, split, scenario and method (2 x 9 x 5)", nrow(k) == 2 * 9 * 5)))
cat("\n")
for (x in checks) cat(if (isTRUE(x[[2]])) "ok  " else "FAIL", x[[1]], "\n")
if (!all(vapply(checks, function(x) isTRUE(x[[2]]), logical(1)))) quit(status = 1)
cat("all synthetic checks passed\n")
