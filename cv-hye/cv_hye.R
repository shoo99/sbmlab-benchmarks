## Replicate CV on the PXD028735 HYE benchmark (3 injections per condition), three engines.
## CV = sd/mean of linear intensities, per protein and condition, only when all 3 replicates have a value.
## Also: MaxQuant raw "Intensity" vs "LFQ intensity", CV by abundance quartile, and whether low CV predicts an accurate ratio.
## Truth log2(B/A): human 0, yeast -1, E. coli +2. Accuracy = |centred log2 ratio - truth| <= 0.5 (centred on the human median).
suppressPackageStartupMessages(library(data.table))
con <- sub("^>", "", sub(" .*", "", grep("^>", readLines("contaminants_maxquant_2.8.1.fasta"), value = TRUE)))
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }
spmq <- function(h) { ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox]); s <- s[!is.na(s)]
  if (length(s) == 1) unname(s) else NA_character_ }
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
spq <- vapply(pg$`Fasta headers`, spmq, "")
cp <- fread("/home/sysoft/fp-hye/work/combined_protein.tsv", sep = "\t", quote = "")
acc_fp <- function(p, ind) c(sub("^sp\\|([^|]+)\\|.*", "\\1", p), unlist(regmatches(ind, gregexpr("(?<=sp\\|)[^|]+", ind, perl = TRUE))))
cp <- cp[!mapply(function(p, i) any(acc_fp(p, i) %in% con), cp$Protein, cp$`Indistinguishable Proteins`)]
dn <- fread("/home/sysoft/blog-pipeline/eval/diann-hye/report.pg_matrix.tsv", sep = "\t")
dn <- dn[!vapply(strsplit(dn$Protein.Group, ";"), function(a) any(a %in% con), TRUE)]
runs <- grep("wiff$", names(dn), value = TRUE)
ord <- sapply(c("A_Sample_Alpha", "A_Sample_Beta", "A_Sample_Gamma", "B_Sample_Alpha", "B_Sample_Beta", "B_Sample_Gamma"), function(k) grep(k, runs, value = TRUE))
E <- list(
  list(name = "MaxQuant LFQ", species = spq, m = as.matrix(pg[, paste("LFQ intensity", c(A, B)), with = FALSE])),
  list(name = "MaxQuant raw Intensity", species = spq, m = as.matrix(pg[, paste("Intensity", c(A, B)), with = FALSE])),
  list(name = "FragPipe MaxLFQ", species = mapply(function(p, e) sp1(paste(p, e)), cp$Protein, cp$`Indistinguishable Proteins`),
       m = as.matrix(cp[, paste(c(A, B), "MaxLFQ Intensity"), with = FALSE])),
  list(name = "DIA-NN pg_matrix", species = vapply(dn$Protein.Names, sp1, ""), m = as.matrix(dn[, ..ord])))
truth <- c(Human = 0, Yeast = -1, Ecoli = 2)
cvf <- function(x) { if (any(is.na(x))) return(NA_real_); sd(x) / mean(x) }
rows <- list(); qrows <- list(); arows <- list()
for (e in E) {
  m <- e$m; m[m == 0] <- NA; s <- e$species; keep <- !is.na(s); m <- m[keep, ]; s <- s[keep]
  cvA <- apply(m[, 1:3], 1, cvf); cvB <- apply(m[, 4:6], 1, cvf)
  for (sp in c("Human", "Yeast", "Ecoli", "All")) {
    i <- if (sp == "All") rep(TRUE, length(s)) else s == sp
    cvs <- c(cvA[i], cvB[i]); cvs <- cvs[!is.na(cvs)]
    rows[[length(rows) + 1]] <- data.table(engine = e$name, species = sp, n_prot = sum(i), n_cv = length(cvs),
      complete_pct = round(100 * length(cvs) / (2 * sum(i)), 1), median_cv = round(100 * median(cvs), 1),
      pct_lt10 = round(100 * mean(cvs < 0.10), 1), pct_lt20 = round(100 * mean(cvs < 0.20), 1))
  }
  ## CV by abundance quartile (condition A, all species, complete rows)
  okA <- !is.na(cvA); mu <- rowMeans(m[okA, 1:3]); q <- cut(mu, quantile(mu, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
  qrows[[length(qrows) + 1]] <- data.table(engine = e$name, quartile = levels(q), median_cv = round(100 * tapply(cvA[okA], q, median), 1),
                                            n = as.integer(table(q)))
  ## does a low CV predict an accurate ratio? complete in both conditions
  ok <- !is.na(cvA) & !is.na(cvB)
  lr <- log2(rowMeans(m[ok, 4:6])) - log2(rowMeans(m[ok, 1:3])); lr <- lr - median(lr[s[ok] == "Human"])
  acc <- abs(lr - truth[s[ok]]) <= 0.5; maxcv <- pmax(cvA[ok], cvB[ok])
  for (sp in c("Human", "Yeast", "Ecoli")) for (band in c("both<20%", "either>=20%")) {
    j <- s[ok] == sp & (if (band == "both<20%") maxcv < 0.20 else maxcv >= 0.20)
    arows[[length(arows) + 1]] <- data.table(engine = e$name, species = sp, cv_band = band, n = sum(j), within_0.5_pct = round(100 * mean(acc[j]), 1))
  }
}
cat("== CV per protein and condition (3 replicates complete)\n"); print(rbindlist(rows))
cat("\n== Median CV by abundance quartile (condition A)\n"); print(rbindlist(qrows))
cat("\n== Ratio accuracy by CV band (max CV of the two conditions)\n"); print(rbindlist(arows))
fwrite(rbindlist(rows), "cv_summary.csv"); fwrite(rbindlist(qrows), "cv_by_quartile.csv"); fwrite(rbindlist(arows), "cv_vs_accuracy.csv")
