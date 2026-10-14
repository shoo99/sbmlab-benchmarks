## Effect of removing MaxQuant-list contaminants on the Figure 1 scores (FragPipe MaxLFQ, DIA-NN pg_matrix).
suppressPackageStartupMessages(library(data.table))
con <- sub("^>", "", sub(" .*", "", grep("^>", readLines("contaminants_maxquant_2.8.1.fasta"), value = TRUE)))
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma"); TR <- c(Human = 0, Yeast = -1, Ecoli = 2)
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }
sc <- function(m, sp, drop) { m[m == 0] <- NA; m <- log2(m); ok <- rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2 & !is.na(sp) & !drop
  fc <- rowMeans(m[, 4:6], na.rm = TRUE) - rowMeans(m[, 1:3], na.rm = TRUE)
  sapply(names(TR), function(s) { i <- ok & sp %in% s; sprintf("%d %.1f%%", sum(i), 100 * mean(abs(fc[i] - TR[s]) <= 0.5)) }) }
cp <- fread("/home/sysoft/fp-hye/work/combined_protein.tsv", sep = "\t", quote = "")
acc <- function(p, ind) c(sub("^sp\\|([^|]+)\\|.*", "\\1", p), unlist(regmatches(ind, gregexpr("(?<=sp\\|)[^|]+", ind, perl = TRUE))))
fpc <- mapply(function(p, i) any(acc(p, i) %in% con), cp$Protein, cp$`Indistinguishable Proteins`)
fsp <- mapply(function(p, e) sp1(paste(p, e)), cp$Protein, cp$`Indistinguishable Proteins`)
fm <- as.matrix(cp[, paste(c(A, B), "MaxLFQ Intensity"), with = FALSE])
cat("FragPipe kept   :", sc(fm, fsp, rep(FALSE, nrow(cp))), "\nFragPipe removed:", sc(fm, fsp, fpc), "\n")
dn <- fread("/home/sysoft/blog-pipeline/eval/diann-hye/report.pg_matrix.tsv", sep = "\t")
runs <- grep("wiff$", names(dn), value = TRUE)
ord <- sapply(c("A_Sample_Alpha", "A_Sample_Beta", "A_Sample_Gamma", "B_Sample_Alpha", "B_Sample_Beta", "B_Sample_Gamma"), function(k) grep(k, runs, value = TRUE))
dnc <- vapply(strsplit(dn$Protein.Group, ";"), function(a) any(a %in% con), TRUE)
dm <- as.matrix(dn[, ..ord]); dsp <- vapply(dn$Protein.Names, sp1, "")
cat("DIA-NN kept     :", sc(dm, dsp, rep(FALSE, nrow(dn))), "\nDIA-NN removed  :", sc(dm, dsp, dnc), "\n")
