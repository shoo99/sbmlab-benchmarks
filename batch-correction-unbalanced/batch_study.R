## Batch correction with unbalanced designs: does correcting first and testing afterwards create false significance?
## (issue #8; Nygaard, Rødland & Hovig 2016, Biostatistics 17:29-39)
##
## Methods. Each one ends in limma lmFit + eBayes; we take the p-value of the group coefficient, BH across genes:
##   ignore batch                         ~ group on the raw data
##   batch in the model                   ~ batch + group
##   ComBat (group protected)             sva::ComBat(mod = ~ group), then ~ group
##   ComBat (batch only)                  sva::ComBat(mod = NULL), then ~ group
##   removeBatchEffect (group protected)  limma::removeBatchEffect(design = ~ group), then ~ group
##
## Part 1, simulation. Two batches of 12 samples, 12 per group in total. Batch 1 holds k samples of group A and
## 12 - k of B, batch 2 the reverse, so k = 6 is balanced and k = 11 nearly confounded. Each gene gets a
## batch-specific location and scale.
##   null:    no group effect, so every call is false
##   planted: 10% of genes change by +-1 log2 between the groups
## Part 2, real data null. bladderbatch: 57 bladder arrays in 5 batches; the 5,000 probes with the highest mean.
## Fake group labels are assigned at random within each batch, with share p of group A in batches 1, 3, 5 and
## 1 - p in batches 2, 4. Conditional on batch the labels carry no biology, so every call is false.
##
## Usage: Rscript batch_study.R > batch_study.out   (about 1-2 min; also writes scores.csv)
suppressPackageStartupMessages({ library(limma); library(sva); library(bladderbatch) })
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
FOLDER <- "batch-correction-unbalanced"
VERSION <- sprintf("limma %s; sva %s", packageVersion("limma"), packageVersion("sva"))
METHODS <- c("ignore batch", "batch in the model", "ComBat (group protected)", "ComBat (batch only)",
             "removeBatchEffect (group protected)")

quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressMessages(expr))); out }
p_group <- function(y, design) eBayes(lmFit(y, design))$p.value[, "groupB"]
p_all_methods <- function(y, batch, group) {
  dg <- model.matrix(~ group)
  db <- model.matrix(~ batch + group)
  list("ignore batch" = p_group(y, dg),
       "batch in the model" = p_group(y, db),
       "ComBat (group protected)" = p_group(quiet(ComBat(y, batch, mod = dg)), dg),
       "ComBat (batch only)" = p_group(quiet(ComBat(y, batch, mod = NULL)), dg),
       "removeBatchEffect (group protected)" = p_group(removeBatchEffect(y, batch = batch, design = dg), dg))
}
## per repetition: share p < 0.05 and p < 0.01 among null genes, BH discoveries, any false BH discovery, power, FDP
rep_stats <- function(p, truth) {
  bh <- p.adjust(p, "BH") < 0.05
  c(p05 = mean(p[!truth] < 0.05), p01 = mean(p[!truth] < 0.01), disc = sum(bh), any_false = as.numeric(any(bh & !truth)),
    power = if (any(truth)) mean(bh[truth]) else NA, fdp = sum(bh & !truth) / max(sum(bh), 1))
}
table_of <- function(res, stat, digits) {
  t(sapply(METHODS, function(m) sapply(names(res), function(k) round(mean(res[[k]][[m]][, stat]), digits))))
}
rows <- list()
add <- function(dataset, tool, setting, subset, metric, value, n) {
  rows[[length(rows) + 1]] <<- scores_rows(FOLDER, dataset, tool, VERSION, setting, subset, metric, value, n)
}

## ---- Part 1: simulation
set.seed(2026)
G <- 1000; NB <- 12; REPS <- 40; K <- c(6, 8, 9, 10, 11)
simulate <- function(k, planted) {
  batch <- factor(rep(1:2, each = NB))
  group <- factor(c(rep("A", k), rep("B", NB - k), rep("A", NB - k), rep("B", k)))
  b <- as.integer(batch)
  loc <- matrix(rnorm(G * 2, 0, 0.6), G)                 # batch location per gene
  sca <- matrix(exp(rnorm(G * 2, 0, 0.25)), G)           # batch scale per gene
  sigma <- exp(rnorm(G, log(0.4), 0.3))                  # residual SD per gene
  truth <- rep(FALSE, G); beta <- numeric(G)
  if (planted) { i <- sample(G, G / 10); truth[i] <- TRUE; beta[i] <- sample(c(-1, 1), length(i), TRUE) }
  y <- rnorm(G, 8, 1.5) + loc[, b] + sca[, b] * sigma * matrix(rnorm(G * 2 * NB), G) + outer(beta, group == "B")
  list(y = y, batch = batch, group = group, truth = truth)
}
cat(sprintf("== Part 1: simulation, %d genes, 2 batches x %d samples (12 per group), %d repetitions per setting\n", G, NB, REPS))
cat("   k = samples of group A in batch 1 (batch 2 holds 12 - k); k = 6 is balanced\n")
for (scenario in c("null", "planted")) {
  res <- setNames(lapply(K, function(k) {
    r <- lapply(seq_len(REPS), function(i) {
      d <- simulate(k, scenario == "planted")
      lapply(p_all_methods(d$y, d$batch, d$group), rep_stats, truth = d$truth)
    })
    setNames(lapply(METHODS, function(m) do.call(rbind, lapply(r, `[[`, m))), METHODS)
  }), paste0("k=", K))
  dataset <- sprintf("simulation (%s): %d genes, 2 batches x %d, %d repetitions", scenario, G, NB, REPS)
  if (scenario == "null") {
    cat("\n-- null (no group effect): share of genes with p < 0.05 (expected 0.05)\n"); print(table_of(res, "p05", 3))
    cat("\n-- null: share with p < 0.01 (expected 0.01)\n"); print(table_of(res, "p01", 3))
    cat("\n-- null: BH-significant genes per repetition (all false; expected about 0)\n"); print(table_of(res, "disc", 1))
    cat("\n-- null: share of repetitions with any BH-significant gene (BH at 5% keeps this <= 0.05 under the null)\n"); print(table_of(res, "any_false", 3))
    for (k in names(res)) for (m in METHODS) {
      s <- res[[k]][[m]]; setting <- sprintf("batch 1 = %s A + %d B", sub("k=", "", k), NB - as.integer(sub("k=", "", k)))
      add(dataset, m, setting, "all genes (null)", "share_p_below_0.05", round(mean(s[, "p05"]), 4), G * REPS)
      add(dataset, m, setting, "all genes (null)", "share_p_below_0.01", round(mean(s[, "p01"]), 4), G * REPS)
      add(dataset, m, setting, "all genes (null)", "bh_discoveries_per_rep", round(mean(s[, "disc"]), 2), REPS)
      add(dataset, m, setting, "all genes (null)", "share_reps_any_bh_discovery", round(mean(s[, "any_false"]), 4), REPS)
    }
  } else {
    cat("\n-- planted (10% of genes change by +-1 log2): power at BH 5%\n"); print(table_of(res, "power", 3))
    cat("\n-- planted: false discovery proportion at BH 5% (target 0.05)\n"); print(table_of(res, "fdp", 3))
    cat("\n-- planted: share of unchanged genes with p < 0.05 (expected 0.05)\n"); print(table_of(res, "p05", 3))
    for (k in names(res)) for (m in METHODS) {
      s <- res[[k]][[m]]; setting <- sprintf("batch 1 = %s A + %d B", sub("k=", "", k), NB - as.integer(sub("k=", "", k)))
      add(dataset, m, setting, "changed genes", "power_bh_0.05", round(mean(s[, "power"]), 4), G / 10 * REPS)
      add(dataset, m, setting, "BH discoveries", "fdp_bh_0.05", round(mean(s[, "fdp"]), 4), REPS)
      add(dataset, m, setting, "unchanged genes", "share_p_below_0.05", round(mean(s[, "p05"]), 4), G * 0.9 * REPS)
    }
  }
}

## ---- Part 2: bladderbatch, fake labels at random within batches
data(bladderdata)
e <- exprs(bladderEset)
batch <- factor(pData(bladderEset)$batch)
e <- e[order(rowMeans(e), decreasing = TRUE)[1:5000], ]
set.seed(2027)
REPS2 <- 30; P <- c(0.5, 0.7, 0.8, 0.9)
cat(sprintf("\n== Part 2: bladderbatch, %d arrays in %d batches (sizes %s), %d probes with the highest mean, %d repetitions\n",
            ncol(e), nlevels(batch), paste(table(batch), collapse = "/"), nrow(e), REPS2))
cat("   p = share of fake group A in batches 1, 3, 5 (1 - p in batches 2, 4); p = 0.5 is balanced\n")
res <- setNames(lapply(P, function(p) {
  r <- lapply(seq_len(REPS2), function(i) {
    group <- character(ncol(e))
    for (b in levels(batch)) {
      idx <- which(batch == b); pa <- if (as.integer(b) %% 2 == 1) p else 1 - p
      group[idx] <- sample(rep(c("A", "B"), c(round(pa * length(idx)), length(idx) - round(pa * length(idx)))))
    }
    lapply(p_all_methods(e, batch, factor(group, levels = c("A", "B"))), rep_stats, truth = rep(FALSE, nrow(e)))
  })
  setNames(lapply(METHODS, function(m) do.call(rbind, lapply(r, `[[`, m))), METHODS)
}), paste0("p=", P))
cat("\n-- share of probes with p < 0.05 (expected 0.05)\n"); print(table_of(res, "p05", 3))
cat("\n-- BH-significant probes per repetition (all false; expected about 0)\n"); print(table_of(res, "disc", 1))
cat("\n-- share of repetitions with any BH-significant probe (BH at 5% keeps this <= 0.05 under the null)\n"); print(table_of(res, "any_false", 3))
dataset <- sprintf("bladderbatch: %d arrays, %d batches, top %d probes, %d repetitions", ncol(e), nlevels(batch), nrow(e), REPS2)
for (p in names(res)) for (m in METHODS) {
  s <- res[[p]][[m]]; setting <- sprintf("share of A in odd batches = %s", sub("p=", "", p))
  add(dataset, m, setting, "all probes (null)", "share_p_below_0.05", round(mean(s[, "p05"]), 4), nrow(e) * REPS2)
  add(dataset, m, setting, "all probes (null)", "bh_discoveries_per_rep", round(mean(s[, "disc"]), 2), REPS2)
  add(dataset, m, setting, "all probes (null)", "share_reps_any_bh_discovery", round(mean(s[, "any_false"]), 4), REPS2)
}

write_scores(do.call(rbind, rows), file.path(script_dir, "scores.csv"))
cat("\n", R.version.string, " | ", VERSION, " | bladderbatch ", format(packageVersion("bladderbatch")), "\n", sep = "")
