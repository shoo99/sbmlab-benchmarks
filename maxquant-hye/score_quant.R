## Score MaxQuant quantification against the known HYE mixing ratios (PXD028735, Van Puyvelde 2022)
## Condition A = Human 65 / Yeast 30 / E. coli 5 ; B = 65 / 15 / 20  ->  log2(B/A): human 0, yeast -1, E. coli +2
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
TXT <- if (length(args)) args[1] else "/home/sysoft/mq-hye/results_mbr_on/txt"   # txt folder of the MBR-on MaxQuant run
pg <- fread(file.path(TXT, "proteinGroups.txt"), sep = "\t", quote = "")
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)

cat("rows in proteinGroups.txt:", nrow(pg), "\n")
stopifnot("Decoy" %in% names(pg))   # MaxQuant 2.8 calls it Decoy (was Reverse)
f_rev  <- pg[["Decoy"]] %in% "+"
f_con  <- pg[["Potential contaminant"]] %in% "+"
f_site <- pg[["Only identified by site"]] %in% "+"
cat(sprintf("Decoy %d · Potential contaminant %d · Only identified by site %d · any of the three %d\n",
            sum(f_rev), sum(f_con), sum(f_site), sum(f_rev | f_con | f_site)))
pg <- pg[!(f_rev | f_con | f_site)]
cat("after the three standard filters:", nrow(pg), "\n")

## species from the FASTA headers of ALL proteins in the group (OX=9606 human, 559292 yeast, 83333 E. coli)
sp_of <- function(h) {
  ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox])
  s <- s[!is.na(s)]
  if (length(s) == 1) s else if (length(s) > 1) "Mixed" else "Unknown"
}
pg[, species := vapply(`Fasta headers`, sp_of, "")]
print(pg[, .N, by = species][order(-N)])

## gene-name semicolons: how many groups carry more than one gene name
cat(sprintf("groups with ';' in Gene names: %d (%.1f%%)\n",
            sum(grepl(";", pg$`Gene names`)), 100 * mean(grepl(";", pg$`Gene names`))))

score <- function(prefix, min_valid = 2) {
  ma <- as.matrix(pg[, paste(prefix, A), with = FALSE]); mb <- as.matrix(pg[, paste(prefix, B), with = FALSE])
  ma[ma == 0] <- NA; mb[mb == 0] <- NA
  la <- log2(ma); lb <- log2(mb)
  ok <- rowSums(!is.na(la)) >= min_valid & rowSums(!is.na(lb)) >= min_valid
  fc <- rowMeans(lb, na.rm = TRUE) - rowMeans(la, na.rm = TRUE)
  ## per-species normalisation is NOT applied: human (constant) median should be 0 if normalisation works
  out <- rbindlist(lapply(names(EXPECT), function(s) {
    i <- ok & pg$species == s; x <- fc[i]; e <- EXPECT[[s]]
    data.table(quant = prefix, species = s, n = sum(i), expected = e,
               median = round(median(x), 3), bias = round(median(x) - e, 3),
               IQR = round(IQR(x), 3), within_0.5 = round(mean(abs(x - e) <= 0.5), 3))
  }))
  ## separability: yeast vs E. coli with a naive cut at +0.5 (halfway between expected human 0 and E. coli 2 is 1)
  out
}
res <- rbind(score("Intensity"), score("iBAQ"), score("LFQ intensity"))
print(res)
fwrite(res, "/home/sysoft/mq-hye/analysis/quant_scores.csv")

## completeness per quant type (all three species, quantified in >=2 of 3 per condition)
for (p in c("Intensity", "iBAQ", "LFQ intensity")) {
  m <- as.matrix(pg[species %in% names(EXPECT), paste(p, c(A, B)), with = FALSE]); m[m == 0] <- NA
  cat(sprintf("%-14s cells with a value: %.1f%% · groups quantified in all 6: %d\n",
              p, 100 * mean(!is.na(m)), sum(rowSums(!is.na(m)) == 6)))
}
