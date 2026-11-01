## The "Complete Minimal Pipeline" from the post, run verbatim on DIA-NN 2.7.0 report.pg_matrix.tsv (PXD028735 SWATH, 3 vs 3),
## then scored against the known ratios: human log2FC 0 (should not be called), yeast -1 (down), E. coli +2 (up).
suppressPackageStartupMessages({library(data.table); library(limma)})
pg <- fread("../../diann-2.7-hye/report.pg_matrix.tsv", sep = "\t")
runs <- grep("Condition_", names(pg), value = TRUE)
runs <- c(grep("Condition_A", runs, value = TRUE), grep("Condition_B", runs, value = TRUE))   # A first = "Control"
data <- as.data.frame(pg[, ..runs]); rownames(data) <- pg$Protein.Group
species <- setNames(ifelse(grepl("_HUMAN", pg$Protein.Names), "Human", ifelse(grepl("_YEAST", pg$Protein.Names), "Yeast",
                    ifelse(grepl("_ECOLI", pg$Protein.Names), "Ecoli", NA))), pg$Protein.Group)

## ---- post code: Complete Minimal Pipeline (unchanged except the input) ----
data[data == 0] <- NA
data_log <- log2(as.matrix(data))
col_med <- apply(data_log, 2, median, na.rm = TRUE)
data_norm <- sweep(data_log, 2, col_med - median(col_med))
group <- factor(c(rep("Control", 3), rep("Treatment", 3)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- eBayes(contrasts.fit(lmFit(data_norm, design),
              makeContrasts(Treatment - Control, levels = design)))
res <- topTable(fit, number = Inf, sort.by = "none")
res$sig <- abs(res$logFC) > 1 & res$adj.P.Val < 0.05
res <- res[!is.na(res$adj.P.Val), ]
cat(sprintf("%d significant (%d up, %d down) out of %d tested\n",
    sum(res$sig), sum(res$sig & res$logFC > 0), sum(res$sig & res$logFC < 0), nrow(res)))
## ---- end post code ----

score <- function(r, label) {
  r$species <- species[rownames(r)]
  sig <- r$sig; fp <- sig & r$species == "Human"
  wrongdir <- sig & ((r$species == "Yeast" & r$logFC > 0) | (r$species == "Ecoli" & r$logFC < 0))
  cat(sprintf("%s: tested %d | called %d | human (false) %d | yeast down %d/%d | E. coli up %d/%d | wrong direction %d | false share %.1f%%\n",
      label, nrow(r), sum(sig), sum(fp), sum(sig & r$species == "Yeast" & r$logFC < 0), sum(r$species == "Yeast"),
      sum(sig & r$species == "Ecoli" & r$logFC > 0), sum(r$species == "Ecoli"), sum(wrongdir), 100 * (sum(fp) + sum(wrongdir)) / max(1, sum(sig))))
}
score(res, "minimal pipeline (no filter)")
cat("Note: yeast truth is -1, so a |logFC| > 1 cutoff sits exactly at the true yeast change.\n")

## Same with the post's Step 1 per-group filter (>= 2 of 3 in each group) before normalising
keep <- sapply(levels(group), function(g) rowSums(!is.na(data_log[, group == g])) >= 2)
dl <- data_log[apply(keep, 1, all), ]
cm <- apply(dl, 2, median, na.rm = TRUE); dn <- sweep(dl, 2, cm - median(cm))
f2 <- eBayes(contrasts.fit(lmFit(dn, design), makeContrasts(Treatment - Control, levels = design)))
r2 <- topTable(f2, number = Inf, sort.by = "none"); r2$sig <- abs(r2$logFC) > 1 & r2$adj.P.Val < 0.05; r2 <- r2[!is.na(r2$adj.P.Val), ]
score(r2, "with per-group filter (2 of 3)")
## and with a smaller fold-change threshold that can see the yeast change
r3 <- r2; r3$sig <- abs(r3$logFC) > 0.5 & r3$adj.P.Val < 0.05
score(r3, "per-group filter, |logFC| > 0.5")
fwrite(cbind(protein = rownames(r2), species = species[rownames(r2)], as.data.table(r2)), "limma_results_filtered.csv")

## Extra check (added after the first run): the post's median normalization on top of DIA-NN's own normalization.
## In this design more proteins go down (yeast) than up (E. coli), so the column medians shift and re-centring moves everything.
cat(sprintf("\nmedian logFC after the post's median normalization: human %.3f, yeast %.3f, E. coli %.3f\n",
    median(r2$logFC[species[rownames(r2)] %in% "Human"]), median(r2$logFC[species[rownames(r2)] %in% "Yeast"]),
    median(r2$logFC[species[rownames(r2)] %in% "Ecoli"])))
f4 <- eBayes(contrasts.fit(lmFit(dl, design), makeContrasts(Treatment - Control, levels = design)))   # DIA-NN normalization only
r4 <- topTable(f4, number = Inf, sort.by = "none"); r4 <- r4[!is.na(r4$adj.P.Val), ]
cat(sprintf("median logFC without re-normalizing (DIA-NN's own normalization): human %.3f, yeast %.3f, E. coli %.3f\n",
    median(r4$logFC[species[rownames(r4)] %in% "Human"]), median(r4$logFC[species[rownames(r4)] %in% "Yeast"]),
    median(r4$logFC[species[rownames(r4)] %in% "Ecoli"])))
r4$sig <- abs(r4$logFC) > 1 & r4$adj.P.Val < 0.05; score(r4, "per-group filter, no re-normalization, |logFC| > 1")
r5 <- r4; r5$sig <- abs(r5$logFC) > 0.5 & r5$adj.P.Val < 0.05; score(r5, "per-group filter, no re-normalization, |logFC| > 0.5")

## adj.P < 0.05 alone (no fold-change cutoff), with and without the extra median normalization
r6 <- r2; r6$sig <- r6$adj.P.Val < 0.05; score(r6, "per-group filter, extra normalization, adj.P only")
r7 <- r4; r7$sig <- r7$adj.P.Val < 0.05; score(r7, "per-group filter, no re-normalization, adj.P only")
