## Writes synthetic droplets with known doublets in the layout of the two downloads (used by test/test_pipeline.sh).
##   GSE96583 batch-2 layout: 4 donors x 3 cell types, control and stimulated lanes of 1,200 droplets each; 10% of
##     droplets are two cells summed. demuxlet-style labels: "doublet" only when the two donors differ (a same-donor
##     doublet is labelled "singlet", as demuxlet would), and 2% of droplets "ambs".
##   Cell Ranger 2 layout (hgmm_12k/filtered_gene_bc_matrices/{hg19,mm10}): 1,000 droplets, half human, half mouse,
##     8% doublets; each droplet also gets 2% ambient RNA of the other species.
## It also writes generator_truth.tsv: dataset, barcode, doublet (1/0), visible (1/0), for the checks.
## Usage: Rscript test/make_synthetic.R <data_dir>
suppressPackageStartupMessages(library(Matrix))
out <- commandArgs(TRUE)[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(3)
G <- 800
type_mu <- sapply(1:3, function(t) rlnorm(G, -1, 1.6))             # three cell types
cell_counts <- function(n, types, mu = type_mu) {
  m <- sapply(seq_len(n), function(i) mu[, types[i]] * rlnorm(1, log(2), 0.3))
  matrix(rnbinom(length(m), mu = m, size = 3), nrow(mu))
}
barcode <- function(n, tag) sprintf("%s%07d-1", tag, seq_len(n))
gt <- list()

## ---- Kang-like lanes
donors <- c("101", "1015", "1016", "1039")
genes <- data.frame(id = sprintf("ENSG%011d", seq_len(G)), sym = sprintf("GENE%d", seq_len(G)))
write.table(genes, gzfile(file.path(out, "GSE96583_batch2.genes.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
tabs <- list()
for (lane in list(c("ctrl", "GSM2560248", "2.1"), c("stim", "GSM2560249", "2.2"))) {
  n <- 1200; dbl <- runif(n) < 0.10
  t1 <- sample(3, n, TRUE); t2 <- sample(3, n, TRUE)
  d1 <- sample(donors, n, TRUE); d2 <- sample(donors, n, TRUE)
  x <- cell_counts(n, t1); x[, dbl] <- x[, dbl] + cell_counts(sum(dbl), t2[dbl])
  bc <- barcode(n, toupper(substr(lane[1], 1, 1)))
  label <- ifelse(dbl & d1 != d2, "doublet", "singlet")
  label[sample(n, round(0.02 * n))] <- "ambs"
  mtx <- file.path(out, sprintf("%s_%s.mtx", lane[2], lane[3]))
  invisible(writeMM(Matrix(x, sparse = TRUE), mtx)); system2("gzip", c("-nf", shQuote(mtx)))
  writeLines(bc, gzfile(file.path(out, sprintf("%s_barcodes.tsv.gz", lane[2]))))
  tabs[[lane[1]]] <- data.frame(tsne1 = 0, tsne2 = 0, ind = d1, stim = lane[1], cluster = t1,
                                cell = c("CD4 T cells", "B cells", "CD14+ Monocytes")[t1], multiplets = label, row.names = bc)
  gt[[lane[1]]] <- data.frame(dataset = paste0("kang_", lane[1]), barcode = bc, doublet = as.integer(dbl),
                              visible = as.integer(dbl & d1 != d2), label = label)
}
write.table(do.call(rbind, unname(tabs)), gzfile(file.path(out, "GSE96583_batch2.total.tsne.df.tsv.gz")), sep = "\t", quote = FALSE)

## ---- hgmm-like: human and mouse genes, Cell Ranger 2 folders per genome
n <- 1000; dbl <- runif(n) < 0.08
sp1 <- sample(c("human", "mouse"), n, TRUE); sp2 <- sample(c("human", "mouse"), n, TRUE)
types <- sample(3, n, TRUE)
expr <- function(sp, which, idx) {                 # counts of species `which` for droplets idx whose cell is `sp`
  x <- cell_counts(length(idx), types[idx])
  sweep(x, 2, sp[idx] == which, `*`)
}
ambient <- function(k, which) matrix(rpois(G * k, rowMeans(type_mu) * 0.04), G)
bc <- barcode(n, "H")
for (g in c("hg19", "mm10")) {
  which <- if (g == "hg19") "human" else "mouse"
  x <- expr(sp1, which, seq_len(n))
  x[, dbl] <- x[, dbl] + expr(sp2, which, which(dbl))
  x <- x + ambient(n, which)
  dir.create(file.path(out, "hgmm_12k", "filtered_gene_bc_matrices", g), recursive = TRUE, showWarnings = FALSE)
  d <- file.path(out, "hgmm_12k", "filtered_gene_bc_matrices", g)
  invisible(writeMM(Matrix(x, sparse = TRUE), file.path(d, "matrix.mtx")))
  write.table(data.frame(sprintf("%s_G%04d", g, seq_len(G)), sprintf("%s_gene%d", g, seq_len(G))), file.path(d, "genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(bc, file.path(d, "barcodes.tsv"))
}
gt$hgmm <- data.frame(dataset = "hgmm_12k", barcode = bc, doublet = as.integer(dbl), visible = as.integer(dbl & sp1 != sp2),
                      label = NA)
write.table(do.call(rbind, gt), file.path(out, "generator_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("synthetic data in %s: Kang-like lanes of 1,200 droplets (10%% doublets), hgmm-like 1,000 droplets (8%%)\n", out))
