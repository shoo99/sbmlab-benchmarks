## Figure 1 scores after centring each engine's ratios on its own human median (removes a global offset).
suppressPackageStartupMessages(library(data.table))
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma"); TR <- c(Human = 0, Yeast = -1, Ecoli = 2)
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }
spmq <- function(h) { ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox]); s <- s[!is.na(s)]; if (length(s) == 1) unname(s) else NA_character_ }
sc <- function(name, m, sp) { m[m == 0] <- NA; m <- log2(m); ok <- rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2 & !is.na(sp)
  fc <- (rowMeans(m[, 4:6], na.rm = TRUE) - rowMeans(m[, 1:3], na.rm = TRUE))[ok]; s <- sp[ok]; off <- median(fc[s == "Human"])
  data.table(engine = name, species = names(TR), n = sapply(names(TR), function(k) sum(s == k)), human_offset = round(off, 3),
             within_uncentred = sapply(names(TR), function(k) round(100 * mean(abs(fc[s == k] - TR[k]) <= 0.5), 1)),
             within_centred = sapply(names(TR), function(k) round(100 * mean(abs(fc[s == k] - off - TR[k]) <= 0.5), 1))) }
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
cp <- fread("/home/sysoft/fp-hye/work/combined_protein.tsv", sep = "\t", quote = "")
dn <- fread("/home/sysoft/blog-pipeline/eval/diann-hye/report.pg_matrix.tsv", sep = "\t")
runs <- grep("wiff$", names(dn), value = TRUE)
ord <- sapply(c("A_Sample_Alpha", "A_Sample_Beta", "A_Sample_Gamma", "B_Sample_Alpha", "B_Sample_Beta", "B_Sample_Gamma"), function(k) grep(k, runs, value = TRUE))
res <- rbind(sc("MaxQuant", as.matrix(pg[, paste("LFQ intensity", c(A, B)), with = FALSE]), vapply(pg$`Fasta headers`, spmq, "")),
             sc("FragPipe", as.matrix(cp[, paste(c(A, B), "MaxLFQ Intensity"), with = FALSE]), mapply(function(p, e) sp1(paste(p, e)), cp$Protein, cp$`Indistinguishable Proteins`)),
             sc("DIA-NN", as.matrix(dn[, ..ord]), vapply(dn$Protein.Names, sp1, "")))
print(res); fwrite(res, "centred_engines.csv")
cm <- fread("/home/sysoft/blog-pipeline/eval/diann-hye/common_3866.csv")
for (col in c("fc_mq", "fc_fp", "fc_dn")) { off <- median(cm[species == "Human", get(col)])
  cat(col, "common-set offset", round(off, 3), "centred within:", sapply(names(TR), function(k) round(100 * mean(abs(cm[species == k, get(col)] - off - TR[k]) <= 0.5), 1)), "\n") }
