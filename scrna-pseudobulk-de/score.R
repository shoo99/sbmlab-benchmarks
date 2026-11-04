## Scores the bundled tables: sim_scores.csv (sim_study.R) and, once the Kang run is done, kang_metrics.csv
## (kang_study.R). Writes scores.csv (simulation rows, then Kang rows) and prints the Kang summary.
##   null      share of genes with p < 0.05 (expected 0.05), BH 5% discoveries per split (all false), and the share of
##             splits with any BH discovery (BH at 5% keeps it <= 0.05)
##   planted   power and FDP of BH 5% calls
## Usage: Rscript score.R > score.out   (seconds; base R and the bundled tables only)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
options(width = 140)
FOLDER <- "scrna-pseudobulk-de"
rows <- list(utils::read.csv(file.path(script_dir, "sim_scores.csv"), stringsAsFactors = FALSE, check.names = FALSE))
kang <- file.path(script_dir, "kang_metrics.csv")
if (file.exists(kang)) {
  k <- utils::read.csv(kang, stringsAsFactors = FALSE, check.names = FALSE)
  version <- utils::read.csv(file.path(script_dir, "kang_versions.csv"), stringsAsFactors = FALSE)$version
  add <- function(dataset, tool, setting, subset, metric, value, n) {
    rows[[length(rows) + 1]] <<- scores_rows(FOLDER, dataset, tool, version, setting, subset, metric, value, n)
  }
  cat("Kang 2018 control PBMC, 8 donors split 4 vs 4\n")
  cat("p05: share of genes with p < 0.05 (null; expected 0.05); false: BH 5% discoveries per split (null, all false);\n")
  cat("any: share of splits with any BH discovery; power / fdp: planted 2-fold changes in 10% of genes, BH 5%\n")
  for (ct in unique(k$cell_type)) {
    d <- k[k$cell_type == ct, ]
    dataset <- sprintf("Kang 2018 control PBMC: %s, %d cells, %d genes", ct, d$cells[1], d$genes[1])
    tab <- list()
    for (m in unique(d$method)) {
      nu <- d[d$method == m & d$scenario == "null", ]
      pl <- d[d$method == m & d$scenario == "planted", ]
      v <- c(p05 = mean(nu$p05), false = mean(nu$false), any = mean(nu$false > 0),
             power = if (nrow(pl)) mean(pl$power) else NA, fdp = if (nrow(pl)) mean(pl$fdp) else NA)
      tab[[m]] <- data.frame(method = m, t(round(v, 3)))
      add(dataset, m, "null: donors split 4 vs 4", "all genes (null)", "share_p_below_0.05", round(v[["p05"]], 4), nrow(nu))
      add(dataset, m, "null: donors split 4 vs 4", "BH 5% discoveries", "false_discoveries_per_rep", round(v[["false"]], 2), nrow(nu))
      add(dataset, m, "null: donors split 4 vs 4", "BH 5% discoveries", "share_reps_any_bh_discovery", round(v[["any"]], 4), nrow(nu))
      if (nrow(pl)) {
        add(dataset, m, "planted 2-fold in 10% (binomial thinning)", "changed genes", "power_bh_0.05", round(v[["power"]], 4), nrow(pl))
        add(dataset, m, "planted 2-fold in 10% (binomial thinning)", "BH 5% discoveries", "fdp_bh_0.05", round(v[["fdp"]], 4), nrow(pl))
      }
    }
    cat(sprintf("\n== %s (%d null splits, %d planted)\n", dataset, length(unique(d$split[d$scenario == "null"])),
                length(unique(d$split[d$scenario == "planted"]))))
    print(do.call(rbind, tab), row.names = FALSE)
  }
} else {
  cat("kang_metrics.csv not present yet (run download.sh and kang_study.R); scores.csv holds the simulation rows only\n")
}
write_scores(do.call(rbind, rows), file.path(script_dir, "scores.csv"))
cat("\n", R.version.string, "\n", sep = "")
