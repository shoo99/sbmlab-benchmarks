suppressPackageStartupMessages({library(data.table)})
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep="\t", quote="")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
data <- as.data.frame(pg[, paste("LFQ intensity", c("A_Alpha","A_Beta","A_Gamma","B_Alpha","B_Beta","B_Gamma")), with=FALSE])
rownames(data) <- sub(";.*", "", pg$`Protein IDs`)
data_raw <- data
set.seed(1)
pdf("plots.pdf")

# data: proteins x samples, raw intensities (e.g. MaxQuant LFQ), 0 = not quantified
data[data == 0] <- NA
data_log <- log2(as.matrix(data))

group <- factor(c(rep("Control", 3), rep("Treatment", 3)))
min_per_group <- 2    # e.g. n - 1 for n = 3
keep <- sapply(levels(group), function(g) rowSums(!is.na(data_log[, group == g])) >= min_per_group)
data_filtered <- data_log[apply(keep, 1, all), ]
cat(sprintf("Kept %d / %d proteins\n", nrow(data_filtered), nrow(data_log)))

# Median normalization: a common default, not the only valid choice
col_medians <- apply(data_filtered, 2, median, na.rm = TRUE)
data_norm <- sweep(data_filtered, 2, col_medians - median(col_medians))

# limma can work with the remaining NAs, so imputation is a choice, not a requirement.
# If you impute, do it on log2 values:
library(imputeLCMD)
data_imputed <- impute.MinProb(data_norm, q = 0.01)   # left-censored (MNAR) assumption

library(limma)
boxplot(data_norm, las = 2, main = "After median normalization")   # are the distributions aligned?
plotMDS(data_norm[complete.cases(data_norm), ], col = as.integer(group))   # do replicates cluster by group?

library(limma)

# Create design matrix
group <- factor(c(rep("Control", 3), rep("Treatment", 3)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Create contrast
contrast <- makeContrasts(Treatment - Control, levels = design)

# Fit model
fit <- lmFit(data_norm, design)
fit2 <- contrasts.fit(fit, contrast)
fit3 <- eBayes(fit2)

# Extract results
results <- topTable(fit3, number = Inf, sort.by = "none")
results$protein_id <- rownames(results)   # e.g. UniProt accessions, if those are your row names

# Add significance annotation
results$significant <- abs(results$logFC) > 1 & results$adj.P.Val < 0.05
results$regulation <- ifelse(!results$significant, "NS",
                      ifelse(results$logFC > 0, "Up", "Down"))

cat(sprintf("Significant: %d (%d up, %d down)\n",
    sum(results$significant),
    sum(results$regulation == "Up"),
    sum(results$regulation == "Down")))

data_norm3 <- cbind(data_norm, data_norm[, 1:3] + matrix(rnorm(nrow(data_norm)*3, 0, 0.2), nrow(data_norm)))

library(limma)
group3 <- factor(c("A","A","A","B","B","B","C","C","C"))
design <- model.matrix(~ 0 + group3); colnames(design) <- levels(group3)
contr <- makeContrasts(BvsA = B - A, CvsA = C - A, CvsB = C - B, levels = design)
fit3g <- eBayes(contrasts.fit(lmFit(data_norm3, design), contr))   # data_norm3: proteins x 9 samples

anyDiff <- topTable(fit3g, number = Inf, sort.by = "none")       # moderated F-test across all contrasts
BvsA    <- topTable(fit3g, coef = "BvsA", number = Inf, sort.by = "none")
summary(decideTests(fit3g, method = "separate", p.value = 0.05, lfc = 1))

library(ggplot2)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#bdc3c7")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(
    title = "Volcano Plot — Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  annotate("text", x = max(results$logFC) * 0.8, y = max(-log10(results$adj.P.Val)) * 0.9,
           label = sprintf("%d Up", sum(results$regulation == "Up")), color = "#e74c3c") +
  annotate("text", x = min(results$logFC) * 0.8, y = max(-log10(results$adj.P.Val)) * 0.9,
           label = sprintf("%d Down", sum(results$regulation == "Down")), color = "#3498db")

# UniProt's own accession pattern (a plain [A-Z][0-9A-Z]{5} also matches gene symbols like COL1A1)
uniprot_re <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$"
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
hs <- pg$`Protein IDs`[grepl("OX=9606", pg$`Fasta headers`)]; hs <- sub(";.*","",hs)
sig_ids <- head(intersect(results$protein_id[order(results$P.Value)], hs), 300)   # TEST HARNESS: human IDs (HYE human proteins are unchanged by design)
if (all(grepl(uniprot_re, sig_ids))) {
  ego <- enrichGO(gene = sub("-[0-9]+$", "", sig_ids),   # drop isoform suffixes
                  OrgDb = org.Hs.eg.db, keyType = "UNIPROT", ont = "BP")
}

# WRONG — thousands of false positives
significant <- results$P.Value < 0.05

# RIGHT — BH-corrected
significant <- results$adj.P.Val < 0.05

# WRONG — large fold changes can be noise
significant <- abs(results$logFC) > 1

# RIGHT — both thresholds
significant <- abs(results$logFC) > 1 & results$adj.P.Val < 0.05

# 1. Load & preprocess
data <- data_raw
data[data == 0] <- NA
data_log <- log2(as.matrix(data))
col_med <- apply(data_log, 2, median, na.rm = TRUE)
data_norm <- sweep(data_log, 2, col_med - median(col_med))

# 2. limma DE
library(limma)
group <- factor(c(rep("Control", 3), rep("Treatment", 3)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- eBayes(contrasts.fit(lmFit(data_norm, design),
              makeContrasts(Treatment - Control, levels = design)))
res <- topTable(fit, number = Inf, sort.by = "none")
res$sig <- abs(res$logFC) > 1 & res$adj.P.Val < 0.05

# 3. Results
res <- res[!is.na(res$adj.P.Val), ]   # drop proteins with too few values to test
cat(sprintf("%d significant (%d up, %d down) out of %d tested\n",
    sum(res$sig),
    sum(res$sig & res$logFC > 0),
    sum(res$sig & res$logFC < 0),
    nrow(res)))

cat("\nOK: results", nrow(results), "sig", sum(results$significant), " | F-test rows", nrow(anyDiff), " | minimal sig", sum(res$sig, na.rm=TRUE), " | uniprot ok", exists("ego"), "\n")
