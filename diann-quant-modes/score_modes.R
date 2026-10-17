## Score DIA-NN protein quantities from different quantification routes on the HYE known ratios (PXD028735, SWATH, 3 vs 3).
## All routes in one folder come from the same .quant files (see build_inputs.py for the sets). Truth log2(B/A): human 0, yeast -1, E. coli +2.
## Rules as in the earlier posts: >= 2 of 3 values per condition, ratio = mean log2(B) - mean log2(A), within +-0.5 of truth
## (uncentred; centred = after subtracting the human median). CV = sd/mean of 3 linear values per condition (3/3 only).
## Usage: Rscript score_modes.R [dir]   (dir = clean/nombr or clean/mbr for the clean rerun; default this folder). Missing routes are skipped.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); if (length(a)) setwd(a[1])
S <- c("A_Alpha", "A_Beta", "A_Gamma", "B_Alpha", "B_Beta", "B_Gamma")
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }
con <- sub("^>", "", sub(" .*", "", grep("^>", readLines("/home/sysoft/mq-hye/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta"), value = TRUE)))
names_map <- unique(fread("pg_quantums_precision.tsv")[, .(Protein.Group, Protein.Names)])
names_map <- unique(rbind(names_map, unique(fread("pg_precsum_quantums_precision.tsv")[, .(Protein.Group, Protein.Names)])))
rd <- function(f) if (file.exists(f)) fread(f) else NULL
dlf <- "fragments_wide.aq_reformat.tsv.protein_intensities.tsv"
dl <- if (file.exists(dlf)) { x <- fread(dlf); setnames(x, "protein", "Protein.Group"); merge(names_map, x, by = "Protein.Group") }
M <- list(
  "QuantUMS high-precision (default pg_matrix)" = rd("pg_quantums_precision.tsv"),
  "QuantUMS high-accuracy (pg_matrix)"          = rd("pg_quantums_accuracy.tsv"),
  "Legacy direct MaxLFQ (pg_matrix)"            = rd("pg_legacy_direct.tsv"),
  "Precursor sum (QuantUMS precursors)"         = rd("pg_precsum_quantums_precision.tsv"),
  "Precursor sum (legacy precursors)"           = rd("pg_precsum_legacy_direct.tsv"),
  "Precursor sum, precursors in all 6 runs (QuantUMS)" = rd("pg_precsum6_quantums_precision.tsv"),
  "Fragment MaxLFQ (iq)"                        = rd("pg_fragment_maxlfq_iq.tsv"),
  "Fragment directLFQ"                          = dl)
M <- M[!vapply(M, is.null, TRUE)]
truth <- c(Human = 0, Yeast = -1, Ecoli = 2)
prep <- function(t) {
  t <- t[!vapply(strsplit(t$Protein.Group, ";"), function(a) any(a %in% con), TRUE)]
  sp <- vapply(t$Protein.Names, sp1, ""); m <- as.matrix(t[, ..S]); m[!is.finite(m) | m <= 0] <- NA
  ok <- !is.na(sp) & rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2
  list(pg = t$Protein.Group[ok], sp = sp[ok], m = m[ok, , drop = FALSE])
}
P <- lapply(M, prep)
## the all-6-run precursor sum only covers proteins with complete precursors by construction, so it is left out of the
## common set (its rows in the "same proteins" table are the common proteins it covers; see n)
common <- Reduce(intersect, lapply(P[!grepl("all 6 runs", names(P))], `[[`, "pg"))
cat("proteins scoreable by every route:", length(common), "\n\n")
cvf <- function(x) if (any(is.na(x))) NA_real_ else sd(x) / mean(x)
score <- function(p, keep = NULL) {
  i <- if (is.null(keep)) seq_along(p$pg) else which(p$pg %in% keep)
  lm <- log2(p$m[i, , drop = FALSE]); sp <- p$sp[i]
  lr <- rowMeans(lm[, 4:6], na.rm = TRUE) - rowMeans(lm[, 1:3], na.rm = TRUE); lc <- lr - median(lr[sp == "Human"])
  cv <- c(apply(p$m[i, 1:3, drop = FALSE], 1, cvf), apply(p$m[i, 4:6, drop = FALSE], 1, cvf))
  rbindlist(lapply(names(truth), function(s) { j <- sp == s
    data.table(species = s, n = sum(j), median_log2 = round(median(lr[j]), 3), iqr = round(IQR(lr[j]), 3),
               within_uncentred = round(100 * mean(abs(lr[j] - truth[s]) <= 0.5), 1),
               within_centred = round(100 * mean(abs(lc[j] - truth[s]) <= 0.5), 1),
               median_cv = round(100 * median(c(cv[seq_along(sp)][j], cv[length(sp) + seq_along(sp)][j]), na.rm = TRUE), 1)) }))
}
all <- rbindlist(lapply(names(P), function(k) cbind(route = k, score(P[[k]]))))
com <- rbindlist(lapply(names(P), function(k) cbind(route = k, score(P[[k]], common))))
options(width = 200)
cat("== Each route, all proteins it can score\n"); print(all)
cat("\n== Same proteins for every route (n =", length(common), ")\n"); print(com)
fwrite(all, "scores_all.csv"); fwrite(com, "scores_common.csv")
