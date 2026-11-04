## Survival tests of one continuous marker, shared by cutpoint_sim.R and tcga_signatures.R (needs survival, maxstat).
## marker_tests(time, event, x) returns
##   p   the p-value of each analysis in ANALYSES
##   hr  for the two splits, the log-rank hazard ratio of high vs low, (O/E high) / (O/E low), folded to >= 1
## cox_p(time, event, x, covariates) is the Wald p of x in a Cox model with optional adjustment covariates.
ANALYSES <- c("Cox, continuous score", "median split", "optimal cut-point, naive p", "optimal cut-point, corrected p")
MINPROP <- 0.1   # the cut-point search runs between the 10th and 90th percentile (survminer's default)

cox_p <- function(time, event, x, covariates = NULL) {
  d <- data.frame(time = time, event = event, x = x)
  if (!is.null(covariates)) d <- cbind(d, covariates)
  # warnings about a diverging adjustment coefficient (a sparse stage level) do not affect the p-value of x
  fit <- tryCatch(suppressWarnings(survival::coxph(survival::Surv(time, event) ~ ., data = d)), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  unname(summary(fit)$coefficients["x", "Pr(>|z|)"])
}

split_test <- function(time, event, high) {
  sd <- survival::survdiff(survival::Surv(time, event) ~ high)
  oe <- sd$obs / sd$exp                    # groups in the order FALSE, TRUE
  hr <- oe[2] / oe[1]
  c(p = stats::pchisq(sd$chisq, 1, lower.tail = FALSE), hr = unname(max(hr, 1 / hr)))
}

marker_tests <- function(time, event, x) {
  med <- split_test(time, event, x > stats::median(x))
  ms <- tryCatch(maxstat::maxstat.test(survival::Surv(time, event) ~ x, data = data.frame(time, event, x),
                                       smethod = "LogRank", pmethod = "Lau92", minprop = MINPROP, maxprop = 1 - MINPROP),
                 error = function(e) NULL)
  opt <- if (is.null(ms)) c(p = NA, hr = NA) else split_test(time, event, x > ms$estimate)
  p <- c(cox_p(time, event, x), med[["p"]], opt[["p"]], if (is.null(ms)) NA else min(1, ms$p.value))
  list(p = stats::setNames(p, ANALYSES),
       hr = c("median split" = med[["hr"]], "optimal cut-point, naive p" = opt[["hr"]]))
}
