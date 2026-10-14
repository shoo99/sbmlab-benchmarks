## t-test vs limma vs DEqMS at n = 3 vs 3, scored against the known HYE ratios (PXD028735).
## MaxQuant 2.8.1 proteinGroups.txt with MBR (same run as the sbmlab MaxQuant benchmark post).
## Truth: log2(B/A) human 0 (unchanged), yeast -1, E. coli +2 (changed).
suppressPackageStartupMessages({library(data.table); library(limma); library(DEqMS)})
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
sp_of <- function(h) { ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox]); s <- s[!is.na(s)]
  if (length(s) == 1) s else "Mixed/Unknown" }
pg[, species := vapply(`Fasta headers`, sp_of, "")]
pg <- pg[species %in% c("Human", "Yeast", "Ecoli")]
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
m <- as.matrix(pg[, paste("LFQ intensity", c(A, B)), with = FALSE]); m[m == 0] <- NA; m <- log2(m)
rownames(m) <- pg$`Protein IDs`
grp <- factor(rep(c("A", "B"), each = 3))
keep <- rowSums(!is.na(m[, 1:3])) == 3 & rowSums(!is.na(m[, 4:6])) == 3   # complete in both groups: every method sees the same proteins
m <- m[keep, ]; info <- pg[keep]
cat("proteins tested:", nrow(m), "\n"); print(table(info$species))

## Welch t-test
tt <- apply(m, 1, function(x) t.test(x[4:6], x[1:3])$p.value)
## limma
design <- model.matrix(~ grp); fit <- eBayes(lmFit(m, design))
lim <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
## DEqMS: variance prior by peptide count (vignette: min razor+unique peptides across samples, +1)
pep <- as.matrix(info[, paste("Razor + unique peptides", c(A, B)), with = FALSE])
fit$count <- apply(pep, 1, min) + 1
deq <- spectraCounteBayes(fit)
dq <- outputResult(deq, coef_col = 2); dq <- dq[rownames(m), ]

lfc <- rowMeans(m[, 4:6]) - rowMeans(m[, 1:3])
single <- info$`Razor + unique peptides` <= 1
res <- data.frame(species = info$species, lfc = lfc, single = single,
                  p_t = tt, adj_t = p.adjust(tt, "BH"), adj_limma = lim$adj.P.Val, adj_deqms = dq$sca.adj.pva)
score <- function(adj, fcmin) {
  call <- adj < 0.05 & abs(res$lfc) >= fcmin
  tp <- sum(call & res$species != "Human"); fp <- sum(call & res$species == "Human")
  pos <- sum(res$species != "Human")
  sprintf("called %4d | true (yeast+E.coli) %4d of %4d (%.0f%%) | false (human) %3d | FDP %.1f%% | single-peptide among calls %d",
          sum(call), tp, pos, 100 * tp / pos, fp, 100 * fp / max(sum(call), 1), sum(call & res$single))
}
for (fcmin in c(0, 0.5)) {
  cat(sprintf("\n== adj.p < 0.05, |log2FC| >= %.1f ==\n", fcmin))
  cat("Welch t   :", score(res$adj_t, fcmin), "\n")
  cat("limma     :", score(res$adj_limma, fcmin), "\n")
  cat("DEqMS     :", score(res$adj_deqms, fcmin), "\n")
}
cat("\nsingle-peptide proteins tested:", sum(res$single), "\n")
fwrite(res, "limma_deqms_hye_results.csv")
