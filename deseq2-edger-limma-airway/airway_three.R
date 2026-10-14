## DESeq2 vs edgeR (QL) vs limma-voom on Bioconductor 'airway' (Himes et al. 2014, GSE52778):
## 4 airway smooth-muscle cell lines, untreated vs dexamethasone. The post's code, as written (~condition),
## and the paired design (~cell + condition) that the experiment actually calls for.
suppressPackageStartupMessages({library(airway); library(DESeq2); library(edgeR); library(limma)})
data(airway)
counts <- assay(airway)
coldata <- as.data.frame(colData(airway))[, c("cell", "dex")]
coldata$condition <- factor(ifelse(coldata$dex == "trt", "treated", "control"), levels = c("control", "treated"))
cat("genes", nrow(counts), "samples", ncol(counts), "\n"); print(table(coldata$cell, coldata$condition))

run <- function(form, label) {
  t0 <- Sys.time()
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = form)
  keep <- rowSums(counts(dds) >= 10) >= 3; dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE); res <- results(dds, contrast = c("condition", "treated", "control"))
  t_de <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dge <- DGEList(counts = counts, group = coldata$condition)
  k <- filterByExpr(dge, design = model.matrix(form, data = coldata)); dge <- dge[k, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM"); design <- model.matrix(form, data = coldata)
  dge <- estimateDisp(dge, design); fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = "conditiontreated"); er <- topTags(qlf, n = Inf, sort.by = "none")$table
  t_er <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  v <- voom(dge, design); lf <- eBayes(lmFit(v, design)); lv <- topTable(lf, coef = "conditiontreated", number = Inf, sort.by = "none")
  t_lv <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sets <- function(fc) list(
    DESeq2 = rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) >= fc)],
    edgeR  = rownames(er)[er$FDR < 0.05 & abs(er$logFC) >= fc],
    voom   = rownames(lv)[lv$adj.P.Val < 0.05 & abs(lv$logFC) >= fc])
  cat(sprintf("\n==== %s ====  genes tested: DESeq2 %d, edgeR/voom %d | seconds: DESeq2 %.1f, edgeR %.1f, voom %.1f\n",
              label, nrow(res), nrow(er), t_de, t_er, t_lv))
  for (fc in c(0, 1)) {
    s <- sets(fc); all3 <- Reduce(intersect, s); any <- Reduce(union, s)
    cat(sprintf("FDR<0.05 & |log2FC|>=%d: DESeq2 %d, edgeR %d, voom %d | in all three %d | in any %d | all-three/any %.0f%%\n",
                fc, length(s$DESeq2), length(s$edgeR), length(s$voom), length(all3), length(any), 100 * length(all3) / length(any)))
    cat(sprintf("   pairwise overlap: D&E %d, D&V %d, E&V %d\n", length(intersect(s$DESeq2, s$edgeR)),
                length(intersect(s$DESeq2, s$voom)), length(intersect(s$edgeR, s$voom))))
  }
  invisible(sets(0))
}
a <- run(~ condition, "as written in the post: ~ condition (pairing ignored)")
b <- run(~ cell + condition, "paired: ~ cell + condition")
cat("\nDESeq2 DEGs gained by modelling the cell line:", length(setdiff(b$DESeq2, a$DESeq2)), " lost:", length(setdiff(a$DESeq2, b$DESeq2)), "\n")
sessionInfo()$running
