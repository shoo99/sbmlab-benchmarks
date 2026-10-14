## MaxQuant LFQ vs FragPipe MaxLFQ (both Orbitrap DDA) vs DIA-NN pg_matrix (TripleTOF SWATH), all with MBR,
## on the proteins all three quantified (matched on leading UniProt accession). Same scoring rules as the other posts.
suppressPackageStartupMessages(library(data.table))
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)
fcs <- function(m, A, B) { m[m == 0] <- NA; m <- log2(m)
  ok <- rowSums(!is.na(m[, A])) >= 2 & rowSums(!is.na(m[, B])) >= 2
  list(ok = ok, fc = rowMeans(m[, B], na.rm = TRUE) - rowMeans(m[, A], na.rm = TRUE)) }
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
mq <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
mq <- mq[!(Decoy %in% "+") & !(`Potential contaminant` %in% "+") & !(`Only identified by site` %in% "+")]
mq[, acc := sub(";.*", "", `Majority protein IDs`)]
m <- as.matrix(mq[, paste("LFQ intensity", c(A, B)), with = FALSE]); colnames(m) <- c(A, B)
r <- fcs(m, A, B); mq[, `:=`(ok = r$ok, fc_mq = r$fc)]
fp <- fread("/home/sysoft/fp-hye/work/combined_protein.tsv", sep = "\t", quote = "")
fp[, acc := `Protein ID`]
m <- as.matrix(fp[, paste(c(A, B), "MaxLFQ Intensity"), with = FALSE]); colnames(m) <- c(A, B)
r <- fcs(m, A, B); fp[, `:=`(ok = r$ok, fc_fp = r$fc)]
dn <- fread("report.pg_matrix.tsv", sep = "\t")
dn[, acc := sub(";.*", "", Protein.Group)]
dn[, species := fifelse(grepl("_HUMAN", Protein.Names) & !grepl("_YEAST|_ECOLI", Protein.Names), "Human",
                fifelse(grepl("_YEAST", Protein.Names) & !grepl("_HUMAN|_ECOLI", Protein.Names), "Yeast",
                fifelse(grepl("_ECOLI", Protein.Names) & !grepl("_HUMAN|_YEAST", Protein.Names), "Ecoli", NA_character_)))]
rc <- grep("Condition_", names(dn), value = TRUE)
m <- as.matrix(dn[, rc, with = FALSE]); DA <- grep("Condition_A", rc, value = TRUE); DB <- grep("Condition_B", rc, value = TRUE)
r <- fcs(m, DA, DB); dn[, `:=`(ok = r$ok, fc_dn = r$fc)]
all3 <- Reduce(function(x, y) merge(x, y, by = "acc"),
               list(dn[ok == TRUE & !is.na(species), .(acc, species, fc_dn)], mq[ok == TRUE, .(acc, fc_mq)], fp[ok == TRUE, .(acc, fc_fp)]))
fwrite(all3, "common_3866.csv")
cat("proteins quantified by all three:", nrow(all3), "\n")
w <- function(f, s) round(mean(abs(f - EXPECT[s]) <= 0.5), 3)
print(all3[, .(n = .N, within_MaxQuant = w(fc_mq, species), within_FragPipe = w(fc_fp, species), within_DIANN = w(fc_dn, species),
               IQR_MaxQuant = round(IQR(fc_mq), 3), IQR_FragPipe = round(IQR(fc_fp), 3), IQR_DIANN = round(IQR(fc_dn), 3)), by = species])
only_dn <- dn[ok == TRUE & !is.na(species) & !(acc %in% mq$acc[mq$ok]) & !(acc %in% fp$acc[fp$ok])]
cat("quantified only by DIA-NN:", nrow(only_dn), "\n")
print(only_dn[, .(n = .N, within = w(fc_dn, species)), by = species])
