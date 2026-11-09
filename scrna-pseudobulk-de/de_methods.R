## Differential expression tests for one cell type, shared by sim_study.R and kang_study.R.
## Each function returns one raw p-value per gene (NA when a test cannot be computed); the callers apply BH.
##
## Cell-level tests treat every cell as a replicate (what Seurat's FindMarkers and Scanpy's rank_genes_groups do
## by default); they run on log-normalised expression, log1p(counts / cell total * 10,000):
##   Wilcoxon   rank-sum test, the normal approximation of wilcox.test(exact = FALSE) (the Seurat default)
##   t-test     Welch t-test per gene (Scanpy's default method 't-test' is a Welch t-test)
##   MAST       hurdle model ~ group + cellular detection rate, LRT on group (only when MAST is installed and
##              DE_MAST=1; it is slow and not an Ubuntu package)
## Pseudobulk tests sum the counts of each donor's cells, so the donor is the replicate:
##   edgeR      QL F-test (glmQLFit / glmQLFTest), TMM
##   DESeq2     Wald test, default size factors (p-values before independent filtering)
##   limma-voom voom + lmFit + eBayes, TMM
## Genes are filtered by the caller, once per cell type, so every method tests the same genes.
suppressPackageStartupMessages({ library(Matrix); library(edgeR); library(limma); library(DESeq2) })

CELL_METHODS <- c("Wilcoxon (cells)", "t-test (cells)")
PB_METHODS <- c("edgeR QL (pseudobulk)", "DESeq2 (pseudobulk)", "limma-voom (pseudobulk)")
use_mast <- function() Sys.getenv("DE_MAST") == "1" && requireNamespace("MAST", quietly = TRUE)
de_methods <- function() c(CELL_METHODS, if (use_mast()) "MAST (cells)", PB_METHODS)

quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressMessages(expr))); out }

log_normalise <- function(counts) {
  lib <- Matrix::colSums(counts)
  as.matrix(log1p(t(t(counts) / lib) * 1e4))
}

cell_wilcox <- function(logn, group) {
  # the normal approximation of wilcox.test(x = group B, y = group A, exact = FALSE, correct = TRUE), vectorised over
  # genes: rank sum with the tie correction of the variance and a continuity correction of 0.5
  b <- group == levels(group)[2]
  n1 <- sum(b); n0 <- sum(!b); N <- n1 + n0
  W <- vapply(seq_len(nrow(logn)), function(i) sum(rank(logn[i, ])[b]), numeric(1)) - n1 * (n1 + 1) / 2
  ties <- vapply(seq_len(nrow(logn)), function(i) { t <- rle(sort(logn[i, ]))$lengths; sum(t^3 - t) }, numeric(1))
  sigma <- sqrt(n1 * n0 / 12 * ((N + 1) - ties / (N * (N - 1))))
  z <- W - n1 * n0 / 2
  z <- (z - sign(z) * 0.5) / sigma
  p <- 2 * stats::pnorm(-abs(z))
  p[sigma == 0] <- NA
  p
}

cell_ttest <- function(logn, group) {
  b <- group == levels(group)[2]
  m1 <- rowMeans(logn[, b, drop = FALSE]); m0 <- rowMeans(logn[, !b, drop = FALSE])
  n1 <- sum(b); n0 <- sum(!b)
  v1 <- rowSums((logn[, b, drop = FALSE] - m1)^2) / (n1 - 1)
  v0 <- rowSums((logn[, !b, drop = FALSE] - m0)^2) / (n0 - 1)
  se2 <- v1 / n1 + v0 / n0
  df <- se2^2 / ((v1 / n1)^2 / (n1 - 1) + (v0 / n0)^2 / (n0 - 1))
  p <- 2 * stats::pt(-abs((m1 - m0) / sqrt(se2)), df)
  p[se2 == 0] <- NA
  p
}

cell_mast <- function(logn, group) {
  cdr <- scale(colMeans(logn > 0))
  sca <- MAST::FromMatrix(logn, data.frame(group = group, cdr = as.vector(cdr)), data.frame(primerid = rownames(logn)),
                          check_sanity = FALSE)
  fit <- quiet(MAST::zlm(~ group + cdr, sca))
  lrt <- quiet(MAST::lrTest(fit, "group"))
  lrt[, "hurdle", "Pr(>Chisq)"][rownames(logn)]
}

pseudobulk <- function(counts, donor) {
  donor <- factor(donor)
  pb <- counts %*% Matrix::sparseMatrix(i = seq_along(donor), j = as.integer(donor), x = 1,
                                        dims = c(length(donor), nlevels(donor)))
  pb <- as.matrix(pb); colnames(pb) <- levels(donor)
  pb
}

pb_edger <- function(pb, group) {
  y <- calcNormFactors(DGEList(pb))
  design <- model.matrix(~ group)
  fit <- glmQLFit(estimateDisp(y, design), design)
  glmQLFTest(fit, coef = 2)$table$PValue
}

pb_deseq2 <- function(pb, group) {
  dds <- quiet(DESeqDataSetFromMatrix(round(pb), data.frame(group = group), ~ group))
  dds <- quiet(DESeq(dds))
  results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)$pvalue
}

pb_voom <- function(pb, group) {
  y <- calcNormFactors(DGEList(pb))
  design <- model.matrix(~ group)
  eBayes(lmFit(voom(y, design), design))$p.value[, 2]
}

## run_methods(counts, donor, group): counts genes x cells (sparse or dense) for one cell type; donor and group per
## cell (group a factor with two levels, constant within donor). Returns a genes x methods matrix of p-values.
run_methods <- function(counts, donor, group) {
  group <- droplevels(factor(group))
  logn <- log_normalise(counts)
  pb <- pseudobulk(counts, donor)
  pb_group <- factor(tapply(as.character(group), factor(donor), `[`, 1)[colnames(pb)], levels = levels(group))
  out <- cbind("Wilcoxon (cells)" = cell_wilcox(logn, group), "t-test (cells)" = cell_ttest(logn, group))
  if (use_mast()) out <- cbind(out, "MAST (cells)" = cell_mast(logn, group))
  out <- cbind(out, "edgeR QL (pseudobulk)" = pb_edger(pb, pb_group), "DESeq2 (pseudobulk)" = pb_deseq2(pb, pb_group),
               "limma-voom (pseudobulk)" = pb_voom(pb, pb_group))
  rownames(out) <- rownames(counts)
  out
}

## per-repetition summaries: share of null genes with p < 0.05, BH 5% discoveries, false discoveries, power, FDP
rep_stats <- function(p, truth) {
  p[is.na(p)] <- 1
  bh <- p.adjust(p, "BH") < 0.05
  c(p05 = mean(p[!truth] < 0.05), disc = sum(bh), false = sum(bh & !truth),
    power = if (any(truth)) mean(bh[truth]) else NA, fdp = sum(bh & !truth) / max(sum(bh), 1))
}
