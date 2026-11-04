## Scores the bundled tables: sim_scores.csv (cutpoint_sim.R) and tcga_pvalues.csv (tcga_signatures.R).
## Writes scores.csv (the simulation rows followed by the TCGA rows) and prints the TCGA summary tables.
##   random sets      share of sets with p < 0.05, observed and permuted survival labels
##   fixed signatures the p-value on the observed labels, and the share of permutations with p < 0.05
## Under permuted labels every p < 0.05 is false, so the expected share is 0.05 for a calibrated analysis.
##
## Usage: Rscript score.R > score.out   (seconds; needs only base R and the bundled tables)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
options(width = 140)
FOLDER <- "random-signatures-tcga"
VERSION <- "survival 3.5.8; maxstat 0.7.25"   # the versions tcga_signatures.R ran with (see tcga_signatures.out)
ANALYSES <- c("Cox, continuous score", "median split", "optimal cut-point, naive p", "optimal cut-point, corrected p",
              "Cox, adjusted for proliferation", "Cox, adjusted for proliferation and stage")
SHORT <- c("cox", "median", "opt_naive", "opt_corrected", "cox_adj_prolif", "cox_adj_prolif_stage")
LABELS <- c(observed = "observed survival", permuted = "permuted survival (null)")
RANDOM <- c("random 50-gene sets", "random 50-gene sets, proliferation-free")
FIXED <- c("proliferation (PAM50, 11 genes)", "top 50 proliferation-correlated genes")

sim <- utils::read.csv(file.path(script_dir, "sim_scores.csv"), stringsAsFactors = FALSE, check.names = FALSE)
p <- utils::read.csv(file.path(script_dir, "tcga_pvalues.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(ANALYSES %in% names(p)))
rows <- list(sim)
add <- function(dataset, tool, setting, subset, metric, value, n) {
  rows[[length(rows) + 1]] <<- scores_rows(FOLDER, dataset, tool, VERSION, setting, subset, metric, value, n)
}
share <- function(x) round(mean(x < 0.05, na.rm = TRUE), 4)

cat("Share of signatures (random sets) or permutations (fixed signatures) with p < 0.05;\n")
cat("p for the fixed signatures on the observed labels. Under permuted labels the expected share is 0.05.\n")
cat("cox = Cox, continuous score; median = median split; opt_naive / opt_corrected = optimal cut-point, naive /\n")
cat("corrected p; cox_adj_* = Cox adjusted for proliferation (and stage). NA: not applicable.\n")
for (cohort in unique(p$cohort)) {
  d <- p[p$cohort == cohort, ]
  dataset <- sprintf("TCGA %s, %s: %d patients, %d events", cohort, d$endpoint[1], d$patients[1], d$events[1])
  cat("\n==", dataset, "\n-- share with p < 0.05\n")
  tab <- list(); ptab <- list()
  for (fam in c(RANDOM, FIXED)) for (lab in names(LABELS)) {
    s <- d[d$family == fam & d$labels == lab, ANALYSES]
    if (fam %in% FIXED && lab == "observed") {
      ptab[[length(ptab) + 1]] <- data.frame(signature = fam, t(setNames(signif(unlist(s[1, ]), 2), SHORT)))
      for (a in ANALYSES) if (!is.na(s[[a]])) add(dataset, a, LABELS[[lab]], fam, "p_value", s[[a]], d$patients[1])
      next
    }
    v <- vapply(ANALYSES, function(a) share(s[[a]]), numeric(1))
    v[is.nan(v)] <- NA
    what <- if (fam %in% FIXED) "permutations" else "sets"
    tab[[length(tab) + 1]] <- data.frame(signatures = sprintf("%s (%d %s)", fam, nrow(s), what), labels = lab,
                                         t(setNames(v, SHORT)))
    for (a in ANALYSES) add(dataset, a, LABELS[[lab]], fam, "share_p_below_0.05", v[[a]], sum(!is.na(s[[a]])))
  }
  print(do.call(rbind, tab), row.names = FALSE)
  cat("-- p-value of the fixed signatures, observed labels\n")
  print(do.call(rbind, ptab), row.names = FALSE)
}
out <- do.call(rbind, rows)
out <- out[!is.na(out$value), ]
write_scores(out, file.path(script_dir, "scores.csv"))
cat("\n", R.version.string, "\n", sep = "")
