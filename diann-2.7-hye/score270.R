## Score DIA-NN protein matrices against the known HYE ratios with the same rules used for MaxQuant and FragPipe:
## >=2 of 3 valid per condition, log2(B/A) from mean log2, within +/-0.5 of truth; species must be unambiguous.
suppressPackageStartupMessages(library(data.table))
D <- "."
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)
sp_of <- function(names) {
  s <- unique(regmatches(names, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", names))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]
  if (length(s) == 1) s else NA_character_
}
score <- function(file, label) {
  x <- fread(file.path(D, file), sep = "\t")
  rc <- grep("Condition_", names(x), value = TRUE)
  A <- grep("Condition_A", rc, value = TRUE); B <- grep("Condition_B", rc, value = TRUE)
  stopifnot(length(A) == 3, length(B) == 3)
  x[, species := vapply(Protein.Names, sp_of, "")]
  m <- as.matrix(x[, c(A, B), with = FALSE]); m[m == 0] <- NA; m <- log2(m)
  ok <- rowSums(!is.na(m[, A])) >= 2 & rowSums(!is.na(m[, B])) >= 2
  fc <- rowMeans(m[, B], na.rm = TRUE) - rowMeans(m[, A], na.rm = TRUE)
  hm <- median(fc[ok & x$species %in% "Human"])
  cat(sprintf("%s: %d protein groups, %d species-ambiguous, %.1f%% cells with a value; human median %.3f\n",
              label, nrow(x), sum(is.na(x$species)), 100 * mean(!is.na(m[!is.na(x$species), ])), hm))
  rbindlist(lapply(names(EXPECT), function(s) {
    i <- ok & x$species %in% s; f <- fc[i]; e <- EXPECT[[s]]
    data.table(set = label, species = s, n = sum(i), median = round(median(f), 3), IQR = round(IQR(f), 3),
               within_0.5 = round(mean(abs(f - e) <= 0.5), 3),
               within_0.5_centred = round(mean(abs(f - hm - e) <= 0.5), 3))
  }))
}
res <- rbind(score("report.pg_matrix.tsv", "pg_matrix (MBR, DIA-NN default filter)"),
             score("pg_mbr.tsv", "MBR pass, Q<=0.01 & PG.Q<=0.01"))
print(res)
fwrite(res, file.path(D, "diann_scores_270.csv"))
