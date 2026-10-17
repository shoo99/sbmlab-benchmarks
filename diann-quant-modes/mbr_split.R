## MBR on vs off, split by protein: for each route, proteins scoreable without MBR (clean/nombr, same route) versus
## proteins that only became scoreable with MBR (clean/mbr). Same scoring rules as score_modes.R (>= 2 of 3 per condition,
## uncentred log2(B/A) within +-0.5 of truth, contaminants and multi-species groups removed).
suppressPackageStartupMessages(library(data.table))
S <- c("A_Alpha", "A_Beta", "A_Gamma", "B_Alpha", "B_Beta", "B_Gamma")
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }
con <- sub("^>", "", sub(" .*", "", grep("^>", readLines("/home/sysoft/mq-hye/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta"), value = TRUE)))
truth <- c(Human = 0, Yeast = -1, Ecoli = 2)
R <- c("QuantUMS high-precision (default pg_matrix)" = "pg_quantums_precision.tsv",
       "QuantUMS high-accuracy (pg_matrix)"          = "pg_quantums_accuracy.tsv",
       "Legacy direct MaxLFQ (pg_matrix)"            = "pg_legacy_direct.tsv",
       "Precursor sum (QuantUMS precursors)"         = "pg_precsum_quantums_precision.tsv",
       "Precursor sum, precursors in all 6 runs (QuantUMS)" = "pg_precsum6_quantums_precision.tsv",
       "Fragment MaxLFQ (iq)"                        = "pg_fragment_maxlfq_iq.tsv",
       "Fragment directLFQ"                          = "fragments_wide.aq_reformat.tsv.protein_intensities.tsv")
load1 <- function(dir, f) {
  p <- file.path(dir, f); if (!file.exists(p)) return(NULL); t <- fread(p)
  if ("protein" %in% names(t)) { setnames(t, "protein", "Protein.Group")
    nm <- unique(fread(file.path(dir, "pg_precsum_quantums_precision.tsv"))[, .(Protein.Group, Protein.Names)]); t <- merge(nm, t, by = "Protein.Group") }
  t <- t[!vapply(strsplit(t$Protein.Group, ";"), function(a) any(a %in% con), TRUE)]
  sp <- vapply(t$Protein.Names, sp1, ""); m <- as.matrix(t[, ..S]); m[!is.finite(m) | m <= 0] <- NA
  ok <- !is.na(sp) & rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2
  lm <- log2(m[ok, , drop = FALSE])
  data.table(pg = t$Protein.Group[ok], sp = sp[ok], lr = rowMeans(lm[, 4:6], na.rm = TRUE) - rowMeans(lm[, 1:3], na.rm = TRUE))
}
summ <- function(d) d[, .(n = .N, median_log2 = round(median(lr), 3), iqr = round(IQR(lr), 3),
                          within = round(100 * mean(abs(lr - truth[sp]) <= 0.5), 1)), by = sp]
out <- rbindlist(lapply(names(R), function(k) {
  off <- load1("clean/nombr", R[[k]]); on <- load1("clean/mbr", R[[k]]); if (is.null(off) || is.null(on)) return(NULL)
  on[, group := ifelse(pg %in% off$pg, "also without MBR", "only with MBR")]
  rbind(cbind(route = k, set = "MBR off", summ(off)),
        cbind(route = k, set = "MBR on: also without MBR", summ(on[group == "also without MBR"])),
        cbind(route = k, set = "MBR on: only with MBR", summ(on[group == "only with MBR"])))
}))
options(width = 200); print(out[order(sp, route)]); fwrite(out, "clean/mbr_split.csv")

## Where does the QuantUMS-vs-legacy gap on MBR-added proteins come from? Take the E. coli proteins that QuantUMS default
## only quantified with MBR, split them by whether the all-6-run precursor sum covers them with MBR on (all precursors
## complete) or not, and score every route on exactly those proteins (routes that lack a protein have a smaller n).
k6 <- "Precursor sum, precursors in all 6 runs (QuantUMS)"; kq <- names(R)[1]
qadd <- setdiff(load1("clean/mbr", R[[kq]])[sp == "Ecoli", pg], load1("clean/nombr", R[[kq]])$pg)
c6 <- load1("clean/mbr", R[[k6]])[sp == "Ecoli", pg]
sets <- list("complete precursors" = intersect(qadd, c6), "incomplete precursors" = setdiff(qadd, c6))
same <- rbindlist(lapply(names(sets), function(g) rbindlist(lapply(names(R), function(k) { x <- load1("clean/mbr", R[[k]])[pg %in% sets[[g]]]
  data.table(subset = g, route = k, n = nrow(x), median_log2 = round(median(x$lr), 3), within = round(100 * mean(abs(x$lr - 2) <= 0.5), 1)) }))))
cat("\n== E. coli proteins QuantUMS default only quantified with MBR (n =", length(qadd), "), split by precursor completeness, every route\n")
print(same); fwrite(same, "clean/mbr_added_by_completeness.csv")
