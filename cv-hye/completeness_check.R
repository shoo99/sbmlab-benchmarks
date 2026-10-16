## Are the protein x condition pairs without all 3 LFQ values the low-abundance ones? (MaxQuant, condition A)
## Abundance proxy that exists for every protein: raw "Intensity" mean over the 3 A runs (zeros as missing).
suppressPackageStartupMessages(library(data.table))
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
A <- c("A_Alpha", "A_Beta", "A_Gamma")
lfq <- as.matrix(pg[, paste("LFQ intensity", A), with = FALSE]); raw <- as.matrix(pg[, paste("Intensity", A), with = FALSE])
raw[raw == 0] <- NA; nl <- rowSums(lfq > 0); ab <- log10(rowMeans(raw, na.rm = TRUE)); ok <- is.finite(ab)
q <- cut(ab[ok], quantile(ab[ok], 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
cat("MaxQuant condition A: share of proteins with all 3 LFQ values, by raw-intensity quartile\n")
print(round(100 * tapply(nl[ok] == 3, q, mean), 1))
cat("median log10 raw intensity: complete =", round(median(ab[ok & nl == 3]), 2), " incomplete =", round(median(ab[ok & nl < 3]), 2), "\n")

## Proteins measured in all 3 A runs (raw Intensity) but without all 3 LFQ values: why, and how noisy were they?
cvf <- function(x) sd(x) / mean(x)
rawc <- rowSums(!is.na(raw)) == 3
inc <- rawc & nl < 3; com <- rawc & nl == 3
pep <- apply(as.matrix(pg[, paste("Razor + unique peptides", A), with = FALSE]), 1, min)
cat("\nraw Intensity in all 3 A runs:", sum(rawc), "| of these, LFQ incomplete:", sum(inc), "| LFQ = 0 in all 3:", sum(rawc & nl == 0), "\n")
cat("median min razor+unique peptides: LFQ-incomplete =", median(pep[inc]), " LFQ-complete =", median(pep[com]), "\n")
cat("median raw-Intensity CV (condition A): LFQ-incomplete =", round(100 * median(apply(raw[inc, ], 1, cvf)), 1),
    "% (n =", sum(inc), ")  LFQ-complete =", round(100 * median(apply(raw[com, ], 1, cvf)), 1), "% (n =", sum(com), ")\n")
