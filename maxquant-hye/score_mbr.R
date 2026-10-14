## Does match-between-runs add correct or wrong quantities? Peptide-level check against the known ratios.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
TXT <- if (length(args)) args[1] else "/home/sysoft/mq-hye/results_mbr_on/txt"   # txt folder of the MBR-on MaxQuant run
pep <- fread(file.path(TXT, "peptides.txt"), sep = "\t", quote = "")
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)

## accession -> species from the three FASTA files
acc2sp <- rbindlist(lapply(c(Human = "human_9606", Yeast = "yeast_559292", Ecoli = "ecoli_83333"), function(f) {
  h <- grep("^>", readLines(Sys.glob(sprintf("/home/sysoft/mq-hye/fasta/uniprot_sp_%s.fasta", f))), value = TRUE)
  data.table(acc = sub("^>[^|]*\\|([^|]+)\\|.*", "\\1", h))
}), idcol = "species")
pep <- pep[!(Decoy %in% "+") & !(`Potential contaminant` %in% "+")]
pep[, species := acc2sp$species[match(`Leading razor protein`, acc2sp$acc)]]
pep <- pep[!is.na(species)]

ia <- as.matrix(pep[, paste("Intensity", A), with = FALSE]); ib <- as.matrix(pep[, paste("Intensity", B), with = FALSE])
ia[ia == 0] <- NA; ib[ib == 0] <- NA
ta <- as.matrix(pep[, paste("Identification type", A), with = FALSE]); tb <- as.matrix(pep[, paste("Identification type", B), with = FALSE])
full <- rowSums(!is.na(ia)) == 3 & rowSums(!is.na(ib)) == 3          # quantified in all 6 runs
n_match <- rowSums(cbind(ta, tb) == "By matching", na.rm = TRUE)
fc <- rowMeans(log2(ib)) - rowMeans(log2(ia))
## centre on the human median so normalisation does not blur the comparison
fc <- fc - median(fc[full & pep$species == "Human"], na.rm = TRUE)

cat("peptides (no decoy/contaminant, species known):", nrow(pep), "\n")
cat("cells identified by matching:", sum(cbind(ta, tb) == "By matching", na.rm = TRUE),
    "of", sum(cbind(ta, tb) %in% c("By MS/MS", "By matching")), "quantified cells\n")
grp <- ifelse(n_match == 0, "all 6 by MS/MS", ifelse(n_match <= 2, "1-2 runs by matching", "3+ runs by matching"))
res <- rbindlist(lapply(names(EXPECT), function(s) rbindlist(lapply(unique(grp), function(g) {
  i <- full & pep$species == s & grp == g; x <- fc[i]; e <- EXPECT[[s]]
  data.table(species = s, group = g, n = sum(i), median = round(median(x), 3),
             IQR = round(IQR(x), 3), within_0.5 = round(mean(abs(x - e) <= 0.5), 3))
}))))
setorder(res, species, group); print(res)
fwrite(res, "/home/sysoft/mq-hye/analysis/mbr_scores.csv")
