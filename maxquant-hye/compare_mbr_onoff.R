## Protein-level: MBR on vs off, same six runs and settings. Scored against log2(B/A) = 0 / -1 / +2.
suppressPackageStartupMessages(library(data.table))
ON  <- "/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt"
OFF <- "/home/sysoft/mq-hye/raw/combined/txt/proteinGroups.txt"
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)

load_pg <- function(f) {
  pg <- fread(f, sep = "\t", quote = "")
  pg <- pg[!(Decoy %in% "+") & !(`Potential contaminant` %in% "+") & !(`Only identified by site` %in% "+")]
  sp <- function(h) {
    ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
    s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox]); s <- s[!is.na(s)]
    if (length(s) == 1) s else NA_character_
  }
  pg[, species := vapply(`Fasta headers`, sp, "")]
  pg[, key := `Majority protein IDs`]
  pg
}
score <- function(pg, prefix) {
  m <- as.matrix(pg[, paste(prefix, c(A, B)), with = FALSE]); m[m == 0] <- NA; m <- log2(m)
  a <- paste(prefix, A); b <- paste(prefix, B)
  ok <- rowSums(!is.na(m[, a])) >= 2 & rowSums(!is.na(m[, b])) >= 2
  fc <- rowMeans(m[, b], na.rm = TRUE) - rowMeans(m[, a], na.rm = TRUE)
  data.table(pid = pg$key, species = pg$species, ok = ok, fc = fc, cells = rowSums(!is.na(m)))
}
on <- load_pg(ON); off <- load_pg(OFF)
cat("protein groups after filters: MBR on", nrow(on), "· off", nrow(off), "\n")
for (p in c("LFQ intensity", "Intensity")) {
  so <- score(on, p); sf <- score(off, p)
  cat(sprintf("\n== %s ==  cells with a value: on %.1f%% · off %.1f%%\n", p,
              100 * sum(so$cells) / (6 * nrow(so)), 100 * sum(sf$cells) / (6 * nrow(sf))))
  res <- rbindlist(lapply(names(EXPECT), function(s) {
    e <- EXPECT[[s]]; a <- so[ok & species %in% s]; b <- sf[ok & species %in% s]
    data.table(species = s, n_on = nrow(a), n_off = nrow(b),
               within_on = round(mean(abs(a$fc - e) <= 0.5), 3), within_off = round(mean(abs(b$fc - e) <= 0.5), 3),
               median_on = round(median(a$fc), 3), median_off = round(median(b$fc), 3),
               IQR_on = round(IQR(a$fc), 3), IQR_off = round(IQR(b$fc), 3))
  }))
  print(res)
  ## proteins that are scoreable only with MBR on: are they accurate?
  only_on <- so[ok & !(pid %in% sf$pid[sf$ok]) & !is.na(species)]
  cat("scoreable only with MBR on:", nrow(only_on), "\n")
  print(only_on[, .(n = .N, within_0.5 = round(mean(abs(fc - EXPECT[species]) <= 0.5), 3)), by = species])
  fwrite(res, sprintf("/home/sysoft/mq-hye/analysis/mbr_onoff_%s.csv", gsub(" ", "_", p)))
}

## same proteins, both runs: does MBR change the accuracy of proteins quantified either way?
so <- score(on, "LFQ intensity"); sf <- score(off, "LFQ intensity")
both <- merge(so[ok & !is.na(species)], sf[ok & !is.na(species)], by = c("pid", "species"), suffixes = c("_on", "_off"))
cat("\n== LFQ, proteins scoreable in BOTH runs ==\n")
print(both[, .(n = .N, within_on = round(mean(abs(fc_on - EXPECT[species]) <= 0.5), 3),
               within_off = round(mean(abs(fc_off - EXPECT[species]) <= 0.5), 3),
               IQR_on = round(IQR(fc_on), 3), IQR_off = round(IQR(fc_off), 3)), by = species])
