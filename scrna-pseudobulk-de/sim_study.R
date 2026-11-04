## Cell-level vs pseudobulk DE when donors differ: simulation with a known answer (issue #6, part 1)
## (Squair et al. 2021, Nat Commun 12:5692; Zimmerman, Espeland & Langefeld 2021, Nat Commun 12:738)
##
## One cell type from 8 donors, 4 per group. Counts are negative binomial around
##   mean = gene mean x donor effect x cell size factor
##   gene mean     log-normal (meanlog -1, sdlog 1.5) per cell, so most genes are sparse as in droplet data
##   donor effect  log-normal per gene and donor, SD 0.25 (biological variation between people)
##   size factor   log-normal per cell, SD 0.3; dispersion 0.2 per gene (cell-to-cell noise)
## Cells per donor are drawn from 100-400, so donors also differ in cell numbers.
##   null     groups of donors do not differ, so every call is false
##   planted  10% of the genes are 2-fold higher (or lower) in group B, the same in every B donor
## donor SD 0 (no donor effect) is also run as a control: then cells are true replicates and cell-level tests are fine.
## Genes detected in fewer than 5% of cells are dropped before testing (the same genes for every method).
##
## Usage: Rscript sim_study.R > sim_study.out   (about 5 min; also writes sim_scores.csv; SIM_REPS sets the repetitions)
## Seeds: set.seed(600 + setting number) before each setting.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
source(file.path(script_dir, "de_methods.R"))
options(width = 140)
FOLDER <- "scrna-pseudobulk-de"
VERSION <- sprintf("edgeR %s; DESeq2 %s; limma %s; R %s (Wilcoxon rank-sum, Welch t)", packageVersion("edgeR"),
                   packageVersion("DESeq2"), packageVersion("limma"), getRversion())
G <- 2000; DONORS <- 8; REPS <- as.integer(Sys.getenv("SIM_REPS", "10")); FC <- 2

simulate <- function(donor_sd, planted) {
  ncell <- sample(100:400, DONORS, replace = TRUE)
  donor <- factor(rep(sprintf("d%d", seq_len(DONORS)), ncell))
  group <- factor(rep(rep(c("A", "B"), each = DONORS / 2), ncell))
  mu <- rlnorm(G, -1, 1.5)
  eff <- matrix(rlnorm(G * DONORS, 0, donor_sd), G)
  de <- rep(1, G)
  if (planted) {
    idx <- sample(G, G / 10)
    de[idx] <- ifelse(runif(length(idx)) < 0.5, FC, 1 / FC)
  }
  sf <- rlnorm(length(donor), 0, 0.3)
  m <- mu * eff[, as.integer(donor)] * outer(de, as.numeric(group == "B"), function(d, b) ifelse(b == 1, d, 1))
  m <- t(t(m) * sf)
  counts <- matrix(rnbinom(length(m), mu = m, size = 1 / 0.2), G)
  rownames(counts) <- sprintf("g%04d", seq_len(G))
  keep <- rowMeans(counts > 0) >= 0.05
  list(counts = Matrix(counts[keep, ], sparse = TRUE), donor = donor, group = group, truth = de[keep] != 1)
}

SETTINGS <- list(list(sd = 0.25, planted = FALSE, name = "null, donor SD 0.25"),
                 list(sd = 0.25, planted = TRUE, name = "planted 2-fold in 10%, donor SD 0.25"),
                 list(sd = 0, planted = FALSE, name = "null, no donor effect"),
                 list(sd = 0, planted = TRUE, name = "planted 2-fold in 10%, no donor effect"))
rows <- list()
add <- function(tool, setting, subset, metric, value, n) {
  rows[[length(rows) + 1]] <<- scores_rows(FOLDER, sprintf("simulation: %d donors (4 vs 4), %d repetitions", DONORS, REPS),
                                          tool, VERSION, setting, subset, metric, value, n)
}
tables <- list()
for (i in seq_along(SETTINGS)) {
  st <- SETTINGS[[i]]
  set.seed(600 + i)
  res <- lapply(seq_len(REPS), function(r) {
    d <- simulate(st$sd, st$planted)
    p <- run_methods(d$counts, d$donor, d$group)
    rbind(t(apply(p, 2, rep_stats, truth = d$truth)), genes = nrow(p))
  })
  methods <- de_methods()
  s <- t(sapply(methods, function(m) colMeans(do.call(rbind, lapply(res, function(x) x[m, ])))))
  genes <- mean(sapply(res, function(x) x["genes", 1]))
  tab <- data.frame(setting = st$name, method = methods, genes = round(genes),
                    p05 = round(s[, "p05"], 3), disc = round(s[, "disc"], 1), false = round(s[, "false"], 1),
                    power = round(s[, "power"], 3), fdp = round(s[, "fdp"], 3), check.names = FALSE)
  tables[[i]] <- tab
  for (m in methods) {
    subset <- if (st$planted) "unchanged genes" else "all genes (null)"
    add(m, st$name, subset, "share_p_below_0.05", round(s[m, "p05"], 4), round(genes * REPS * if (st$planted) 0.9 else 1))
    add(m, st$name, "BH 5% discoveries", "false_discoveries_per_rep", round(s[m, "false"], 2), REPS)
    if (st$planted) {
      add(m, st$name, "changed genes", "power_bh_0.05", round(s[m, "power"], 4), REPS)
      add(m, st$name, "BH 5% discoveries", "fdp_bh_0.05", round(s[m, "fdp"], 4), REPS)
    }
  }
}
cat(sprintf("== %d donors (4 vs 4), 100-400 cells each, %d genes before filtering, %d repetitions per setting\n", DONORS, G, REPS))
cat("   p05: share of unchanged genes with p < 0.05 (expected 0.05); disc / false: BH 5% discoveries / false ones per\n")
cat("   repetition; power and fdp: share of changed genes found, false share of the discoveries (BH 5%)\n\n")
print(do.call(rbind, tables), row.names = FALSE)
write_scores(do.call(rbind, rows), file.path(script_dir, "sim_scores.csv"))
cat("\n", R.version.string, " | ", VERSION, "\n", sep = "")
