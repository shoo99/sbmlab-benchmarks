## How much does choosing the "optimal" cut-point inflate a survival p-value? (issue #11, part 1: simulation)
##
## A continuous marker (for example a signature score) is turned into a survival p-value in four ways:
##   Cox, continuous score            coxph(Surv ~ score), Wald p
##   median split                     high vs low at the median, log-rank p
##   optimal cut-point, naive p       the cut-point between the 10th and 90th percentile with the largest standardised
##                                    log-rank statistic (maxstat; what survminer::surv_cutpoint picks), then the
##                                    log-rank p of that split as if the cut-point had been fixed in advance
##   optimal cut-point, corrected p   the same search, with the p-value of the maximally selected statistic
##                                    (maxstat, pmethod = "Lau92": Lausen & Schumacher 1992)
## Scenarios. Exponential survival (rate 0.1 per year), uniform censoring over 0-10 years (about 37% events; 42% in the threshold setting).
##   null       the marker is independent of survival, so every p < 0.05 is false (expected share 0.05)
##   linear     log hazard ratio rises linearly with the marker: HR 1.25 per SD
##   threshold  the top 25% of the marker have HR 2, the rest HR 1
## Besides p-values, the hazard ratio of the high vs the low group (folded so that it is >= 1) shows how much the
## optimal split also inflates the effect size under the null (truth: 1).
##
## Usage: Rscript cutpoint_sim.R > cutpoint_sim.out   (under a minute; also writes sim_scores.csv)
## Seeds: set.seed(2026 + scenario number) before each scenario, so scenarios do not depend on each other.
suppressPackageStartupMessages({ library(survival); library(maxstat) })
options(width = 120)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
source(file.path(script_dir, "survival_tests.R"))
FOLDER <- "random-signatures-tcga"
VERSION <- sprintf("survival %s; maxstat %s", packageVersion("survival"), packageVersion("maxstat"))
REPS <- 1000

simulate <- function(n, scenario) {
  x <- rnorm(n)
  lp <- switch(scenario, null = 0, linear = log(1.25) * x, threshold = log(2) * (x > qnorm(0.75)))
  tt <- rexp(n, 0.1 * exp(lp))
  cc <- runif(n, 0, 10)
  list(time = pmin(tt, cc), event = as.integer(tt <= cc), x = x)
}
SCENARIOS <- list(list(n = 100, s = "null"), list(n = 300, s = "null"), list(n = 1000, s = "null"),
                  list(n = 300, s = "linear"), list(n = 300, s = "threshold"))

rows <- list()
add <- function(dataset, tool, setting, subset, metric, value, n) {
  rows[[length(rows) + 1]] <<- scores_rows(FOLDER, dataset, tool, VERSION, setting, subset, metric, value, n)
}
summary_rows <- list()
for (i in seq_along(SCENARIOS)) {
  sc <- SCENARIOS[[i]]
  set.seed(2026 + i)
  res <- t(vapply(seq_len(REPS), function(r) {
    d <- simulate(sc$n, sc$s)
    a <- marker_tests(d$time, d$event, d$x)
    c(a$p, events = mean(d$event), hr_median = a$hr[["median split"]], hr_optimal = a$hr[["optimal cut-point, naive p"]])
  }, numeric(length(ANALYSES) + 3)))
  share <- colMeans(res[, ANALYSES] < 0.05)
  setting <- switch(sc$s, null = "null (no association)", linear = "linear effect, HR 1.25 per SD",
                    threshold = "threshold effect, top 25% HR 2")
  dataset <- sprintf("simulation: n = %d, %d repetitions", sc$n, REPS)
  metric <- if (sc$s == "null") "share_p_below_0.05" else "power_p_below_0.05"
  for (a in ANALYSES) add(dataset, a, setting, "simulated marker", metric, round(share[[a]], 4), REPS)
  if (sc$s == "null") {
    add(dataset, "median split", setting, "simulated marker", "median_folded_hr", round(median(res[, "hr_median"]), 3), REPS)
    add(dataset, "optimal cut-point, naive p", setting, "simulated marker", "median_folded_hr",
        round(median(res[, "hr_optimal"]), 3), REPS)
  }
  summary_rows[[i]] <- data.frame(n = sc$n, scenario = sc$s, events = round(mean(res[, "events"]), 2),
                                  t(round(share, 3)), hr_median = round(median(res[, "hr_median"]), 2),
                                  hr_optimal = round(median(res[, "hr_optimal"]), 2), check.names = FALSE)
}
tab <- do.call(rbind, summary_rows)
names(tab)[match(ANALYSES, names(tab))] <- c("cox", "median", "opt_naive", "opt_corrected")
cat(sprintf("== Share of %d repetitions with p < 0.05 (null: false positives, expected 0.05; otherwise power)\n", REPS))
cat("   cox = Cox, continuous score; median = median split; opt_naive / opt_corrected = optimal cut-point, naive /\n")
cat("   corrected p; hr_*: median hazard ratio of high vs low, folded to >= 1 (truth under the null: 1)\n\n")
print(tab, row.names = FALSE)

write_scores(do.call(rbind, rows), file.path(script_dir, "sim_scores.csv"))
cat("\n", R.version.string, " | ", VERSION, "\n", sep = "")
