## Runs the R doublet detectors on one prepared data set (from prepare.R) and writes scores_r.tsv next to it:
## barcode, tool, score (higher = more doublet-like), call (TRUE / FALSE; NA when the tool makes no call).
##   scDblFinder (default rate)    scDblFinder() with its own rate (1% per 1,000 cells)
##   scDblFinder (expected rate)   dbr = the expected rate in info.tsv
##   scds cxds / bcds / hybrid     co-expression, artificial-doublet classifier, and their combination (scores only)
##   DoubletFinder                 Seurat log-normalisation, 2,000 variable genes, PCs 1-10, pN 0.25, pK from the
##                                 parameter sweep (paramSweep + find.pK, as its authors recommend), nExp = expected
##                                 rate x droplets; score = pANN
## Usage: Rscript run_r_tools.R <work_dir>/<dataset>     env: DF_SWEEP=0 uses pK = 0.09 instead of the sweep
## Seed: set.seed(7) before each tool.
suppressPackageStartupMessages({
  library(Matrix); library(SingleCellExperiment); library(scDblFinder); library(scds)
  library(Seurat); library(DoubletFinder)
})
d <- commandArgs(TRUE)[1]
counts <- as(readMM(gzfile(file.path(d, "counts.mtx.gz"))), "CsparseMatrix")
dimnames(counts) <- list(readLines(file.path(d, "genes.tsv")), readLines(file.path(d, "barcodes.tsv")))
info <- utils::read.delim(file.path(d, "info.tsv"))
rate <- info$expected_rate; n_exp <- round(rate * ncol(counts))
options(pbapply.pb = "none")
if (requireNamespace("pbapply", quietly = TRUE)) invisible(pbapply::pboptions(type = "none"))
quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }
res <- list()
add <- function(tool, score, call = NA) res[[tool]] <<- data.frame(barcode = colnames(counts), tool = tool, score = score, call = call)
timing <- list()
timed <- function(tool, expr) { t0 <- proc.time()[[3]]; force(expr); timing[[tool]] <<- proc.time()[[3]] - t0 }

sce <- SingleCellExperiment(list(counts = counts))
timed("scDblFinder (default rate)", {
  set.seed(7); s <- quiet(scDblFinder(sce))
  add("scDblFinder (default rate)", s$scDblFinder.score, s$scDblFinder.class == "doublet")
})
timed("scDblFinder (expected rate)", {
  set.seed(7); s <- quiet(scDblFinder(sce, dbr = rate))
  add("scDblFinder (expected rate)", s$scDblFinder.score, s$scDblFinder.class == "doublet")
})
timed("scds", {
  set.seed(7); s <- quiet(cxds_bcds_hybrid(sce))     # also runs cxds and bcds and keeps their scores
  add("scds cxds", s$cxds_score); add("scds bcds", s$bcds_score); add("scds hybrid", s$hybrid_score)
})
timed("DoubletFinder", {
  set.seed(7)
  seu <- quiet(CreateSeuratObject(counts))
  seu <- quiet(NormalizeData(seu, verbose = FALSE))
  seu <- quiet(ScaleData(FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE), verbose = FALSE))
  seu <- quiet(RunPCA(seu, npcs = 20, verbose = FALSE))
  pk <<- 0.09
  if (Sys.getenv("DF_SWEEP", "1") != "0") {
    sw <- quiet(summarizeSweep(paramSweep(seu, PCs = 1:10, sct = FALSE), GT = FALSE))
    bc <- quiet(find.pK(sw))
    pk <<- as.numeric(as.character(bc$pK[which.max(bc$BCmetric)]))
  }
  seu <- quiet(doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = pk, nExp = n_exp, sct = FALSE))
  md <- seu@meta.data
  add("DoubletFinder", md[[grep("^pANN_", names(md), value = TRUE)[1]]],
      md[[grep("^DF.classifications_", names(md), value = TRUE)[1]]] == "Doublet")
  cat(sprintf("DoubletFinder pK = %s\n", pk))
})
out <- do.call(rbind, res)
out$score <- signif(out$score, 5)
utils::write.table(out, file.path(d, "scores_r.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
v <- c(scDblFinder = paste("scDblFinder", packageVersion("scDblFinder")), scds = paste("scds", packageVersion("scds")),
       DoubletFinder = sprintf("DoubletFinder %s; Seurat %s", packageVersion("DoubletFinder"), packageVersion("Seurat")))
tools <- unique(out$tool)
utils::write.table(data.frame(tool = tools, version = v[sub(" .*", "", tools)],
                              setting = ifelse(tools == "DoubletFinder", sprintf("pK %s, nExp %d", pk, n_exp), "")),
                   file.path(d, "versions_r.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d droplets, expected rate %.3f; seconds: %s\n", basename(d), ncol(counts), rate,
            paste(sprintf("%s %.0f", names(timing), unlist(timing)), collapse = ", ")))
cat(R.version.string, " | scDblFinder ", format(packageVersion("scDblFinder")), " | scds ", format(packageVersion("scds")),
    " | DoubletFinder ", format(packageVersion("DoubletFinder")), " | Seurat ", format(packageVersion("Seurat")), "\n", sep = "")
