## Is the lower accuracy of matched peptide values due to matching or to low abundance?
## MaxQuant 2.8.1 peptides.txt (MBR on), peptides quantified in all six runs, ratios centred on the human median.
## Within each species, peptides are binned by mean log2 intensity (quartiles of that species) and the
## fraction within ±0.5 of truth is compared between "all 6 by MS/MS" and "any matched" in each bin.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
TXT <- if (length(args)) args[1] else "/home/sysoft/mq-hye/results_mbr_on/txt"
FASTA <- "/home/sysoft/mq-hye/fasta/uniprot_sp_%s.fasta"
pep <- fread(file.path(TXT, "peptides.txt"), sep = "\t", quote = "")
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma"); EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)
acc2sp <- rbindlist(lapply(c(Human = "human_9606", Yeast = "yeast_559292", Ecoli = "ecoli_83333"), function(f) {
  h <- grep("^>", readLines(Sys.glob(sprintf(FASTA, f))), value = TRUE)
  data.table(acc = sub("^>[^|]*\\|([^|]+)\\|.*", "\\1", h)) }), idcol = "species")
pep <- pep[!(Decoy %in% "+") & !(`Potential contaminant` %in% "+")]
pep[, species := acc2sp$species[match(`Leading razor protein`, acc2sp$acc)]]; pep <- pep[!is.na(species)]
ia <- as.matrix(pep[, paste("Intensity", A), with = FALSE]); ib <- as.matrix(pep[, paste("Intensity", B), with = FALSE])
ia[ia == 0] <- NA; ib[ib == 0] <- NA
ty <- as.matrix(pep[, paste("Identification type", c(A, B)), with = FALSE])
full <- rowSums(!is.na(ia)) == 3 & rowSums(!is.na(ib)) == 3
nm <- rowSums(ty == "By matching", na.rm = TRUE)
fc <- rowMeans(log2(ib)) - rowMeans(log2(ia)); fc <- fc - median(fc[full & pep$species == "Human"], na.rm = TRUE)
inten <- rowMeans(log2(cbind(ia, ib)))
d <- data.table(species = pep$species, fc, inten, matched = nm > 0)[full]
d[, bin := cut(inten, quantile(inten, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4)), by = species]
d[, ok := abs(fc - EXPECT[species]) <= 0.5]
res <- d[, .(n = .N, within_0.5 = round(mean(ok), 3)), by = .(species, bin, matched)]
setorder(res, species, bin, matched); print(res)
cat("\nShare of matched peptides per intensity quartile:\n")
print(d[, .(matched_share = round(mean(matched), 3)), by = .(species, bin)][order(species, bin)])
fwrite(res, "mbr_intensity_matched.csv")
