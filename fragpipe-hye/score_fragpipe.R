## Score FragPipe (MSFragger + IonQuant) combined_protein.tsv against the known HYE ratios,
## with the same rules used for MaxQuant: >=2 of 3 valid per condition, log2(B/A) from mean log2, within +/-0.5.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
WORK <- if (length(args)) args[1] else "/home/sysoft/fp-hye/work"
cp <- fread(file.path(WORK, "combined_protein.tsv"), sep = "\t", quote = "")
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)

cat("rows in combined_protein.tsv:", nrow(cp), "\n")
con <- grepl("contam_", cp$Protein)
cat("contaminant rows:", sum(con), "\n")
cp <- cp[!con]
## species: protein plus any indistinguishable proteins must agree
ox <- function(p, extra) {
  s <- unique(regmatches(paste(p, extra), gregexpr("_(HUMAN|YEAST|ECOLI)\\b", paste(p, extra)))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]
  if (length(s) == 1) s else NA_character_
}
cp[, species := mapply(ox, Protein, `Indistinguishable Proteins`)]
print(cp[, .N, by = species])

score <- function(suffix) {
  cols <- paste(c(A, B), suffix)
  stopifnot(all(cols %in% names(cp)))
  m <- as.matrix(cp[, ..cols]); colnames(m) <- c(A, B); m[m == 0] <- NA; m <- log2(m)
  ok <- rowSums(!is.na(m[, A])) >= 2 & rowSums(!is.na(m[, B])) >= 2
  fc <- rowMeans(m[, B], na.rm = TRUE) - rowMeans(m[, A], na.rm = TRUE)
  rbindlist(lapply(names(EXPECT), function(s) {
    i <- ok & cp$species %in% s; x <- fc[i]; e <- EXPECT[[s]]
    data.table(quant = suffix, species = s, n = sum(i), median = round(median(x), 3),
               IQR = round(IQR(x), 3), within_0.5 = round(mean(abs(x - e) <= 0.5), 3))
  }))
}
res <- rbind(score("Intensity"), score("MaxLFQ Intensity"))
print(res)
fwrite(res, file.path(WORK, "..", "fragpipe_scores.csv"))
for (s in c("Intensity", "MaxLFQ Intensity")) {
  m <- as.matrix(cp[!is.na(species), paste(c(A, B), s), with = FALSE]); m[m == 0] <- NA
  cat(sprintf("%-17s cells with a value: %.1f%%\n", s, 100 * mean(!is.na(m))))
}
