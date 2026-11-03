## Microbiome differential abundance: do the methods agree, and which control false positives? (issue #18)
## Nearing et al. (2022, Nat Commun 13:342) found the methods disagree widely on real data. Here the answer is known.
##
## Template: MicrobiomeStat's `smokers` data (Charlson et al. 2010): left-throat 16S OTU counts, one sample per subject,
## samples with >= 1,000 reads, OTUs present in >= 10% of them. Each repetition draws 40 samples and splits them
## 20 / 20 at random, so the groups differ only by what is planted:
##   null      nothing changed: every call is false
##   planted   15 OTUs (present in >= 50% of the template) changed 4-fold up or down in group B, on the absolute scale
##   bloom     as planted, and the most abundant OTU also 8-fold up in group B. The bloom pushes down the relative
##             abundance of every other OTU, though their absolute abundance is unchanged (truth 0).
## Changes are applied to counts, and then every sample, in both groups, is redrawn to its own depth with a
## multinomial. So the data stay compositional, and the truth is absolute.
## Before testing, OTUs present in fewer than 4 of the 40 samples are dropped, for every method alike.
##
## Methods (call = adjusted p < 0.05 unless noted; direction from each method's own effect estimate):
##   ALDEx2           aldex(test = "t", mc.samples = 128, denom = "all"), Welch t, BH (we.eBH)
##   ANCOM-BC2        ancombc2(fix_formula = "group"), q < 0.05 (diff)
##   ANCOM-BC2 + SS   the same run, q < 0.05 and passing ANCOM-BC2's pseudo-count sensitivity check (diff_robust)
##   MaAsLin2         TSS + LOG + LM, BH (q < 0.05, not MaAsLin2's default 0.25)
##   LinDA            linda(feature.dat.type = "count"), padj < 0.05
##   DESeq2           poscounts size factors, Wald, padj < 0.05
##   Wilcoxon CLR     centred log-ratio with pseudo-count 0.5, Wilcoxon rank-sum, BH
##   Wilcoxon TSS     relative abundance, Wilcoxon rank-sum, BH
##
## Usage: OPENBLAS_NUM_THREADS=1 Rscript da_study.R > da_study.out   (env: DA_CORES, default 4; ~15 min on 4 cores)
## Needs the conda environment in this folder (conda-linux-64.lock); writes scores.csv next to the script.
suppressPackageStartupMessages({
  library(MicrobiomeStat); library(ALDEx2); library(ANCOMBC); library(Maaslin2); library(DESeq2); library(parallel)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
FOLDER <- "microbiome-da-methods"
CORES <- as.integer(Sys.getenv("DA_CORES", "4"))
REPS <- 20; N_PER_GROUP <- 20; N_CHANGED <- 15; FOLD <- 4; BLOOM <- 8; MIN_PRESENT <- 4; ELIGIBLE_PREV <- 0.50
quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }
versions <- c(ALDEx2 = "ALDEx2", `ANCOM-BC2` = "ANCOMBC", `ANCOM-BC2 + SS` = "ANCOMBC", MaAsLin2 = "Maaslin2",
              LinDA = "MicrobiomeStat", DESeq2 = "DESeq2", `Wilcoxon CLR` = "stats", `Wilcoxon TSS` = "stats")
METHODS <- names(versions)

## ---- template
data(smokers)
meta0 <- smokers$meta
otu <- smokers$otu[, meta0$AIRWAYSITE == "Throat" & meta0$SIDEOFBODY == "Left"]
otu <- otu[, colSums(otu) >= 1000]
otu <- otu[rowMeans(otu > 0) >= 0.10, ]
rownames(otu) <- paste0("otu", seq_len(nrow(otu)))
eligible <- rownames(otu)[rowMeans(otu > 0) >= ELIGIBLE_PREV]
bloom_otu <- rownames(otu)[which.max(rowSums(otu))]
eligible <- setdiff(eligible, bloom_otu)

## ---- one data set with a known answer
simulate <- function(scenario, rep) {
  set.seed(1000 * match(scenario, c("null", "planted", "bloom")) + rep)
  cnt <- otu[, sample(ncol(otu), 2 * N_PER_GROUP)]
  grp <- factor(rep(c("A", "B"), each = N_PER_GROUP))
  truth <- setNames(rep(0, nrow(cnt)), rownames(cnt))
  if (scenario != "null") {
    ch <- sample(eligible, N_CHANGED)
    truth[ch] <- sample(c(-1, 1), N_CHANGED, TRUE) * log2(FOLD)
  }
  if (scenario == "bloom") truth[bloom_otu] <- log2(BLOOM)
  abs_b <- cnt[, grp == "B"] * 2^truth                                   # absolute change in group B
  mod <- cbind(cnt[, grp == "A"], abs_b)
  depth <- colSums(cnt)
  out <- sapply(seq_len(ncol(mod)), function(j) rmultinom(1, depth[j], mod[, j] / sum(mod[, j])))
  dimnames(out) <- dimnames(cnt)
  keep <- rowSums(out > 0) >= MIN_PRESENT
  list(cnt = out[keep, ], grp = grp, truth = truth[keep])
}

## ---- methods: each returns a data.frame(taxon, p, est, call)
run_methods <- function(cnt, grp, seed) {
  meta <- data.frame(group = grp, sample_id = colnames(cnt), row.names = colnames(cnt))  # 2 columns: see README
  tx <- rownames(cnt)
  set.seed(seed)
  a <- quiet(aldex(as.data.frame(cnt), as.character(grp), mc.samples = 128, test = "t", effect = TRUE, denom = "all", verbose = FALSE))
  b <- quiet(ancombc2(data = cnt, taxa_are_rows = TRUE, meta_data = meta, fix_formula = "group", p_adj_method = "BH",
                      prv_cut = 0, lib_cut = 0, struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, verbose = FALSE, n_cl = 1))$res
  od <- tempfile(); on.exit(unlink(od, recursive = TRUE))
  mz <- quiet(Maaslin2(input_data = as.data.frame(t(cnt)), input_metadata = meta, output = od, min_abundance = 0,
                       min_prevalence = 0, normalization = "TSS", transform = "LOG", analysis_method = "LM",
                       fixed_effects = "group", standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1))$results
  ld <- quiet(linda(feature.dat = cnt, meta.dat = meta, formula = "~group", feature.dat.type = "count", prev.filter = 0,
                    alpha = 0.05, verbose = FALSE))$output$groupB
  dds <- DESeqDataSetFromMatrix(round(cnt), meta, ~group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  ds <- as.data.frame(results(quiet(DESeq(dds, quiet = TRUE)), name = "group_B_vs_A"))
  wil <- function(x) {
    p <- apply(x, 1, function(v) wilcox.test(v[grp == "B"], v[grp == "A"], exact = FALSE)$p.value)
    data.frame(taxon = tx, p = p, est = rowMeans(x[, grp == "B"]) - rowMeans(x[, grp == "A"]), call = p.adjust(p, "BH") < 0.05)
  }
  clr <- apply(log(cnt + 0.5), 2, function(v) v - mean(v))
  tss <- sweep(cnt, 2, colSums(cnt), "/")
  mz <- mz[match(tx, mz$feature), ]; ld <- ld[tx, ]; ds <- ds[tx, ]; a <- a[tx, ]; b <- b[match(tx, b$taxon), ]
  list(
    `ALDEx2` = data.frame(taxon = tx, p = a$we.ep, est = a$diff.btw, call = a$we.eBH < 0.05),
    `ANCOM-BC2` = data.frame(taxon = tx, p = b$p_groupB, est = b$lfc_groupB, call = b$diff_groupB %in% TRUE),
    `ANCOM-BC2 + SS` = data.frame(taxon = tx, p = b$p_groupB, est = b$lfc_groupB, call = b$diff_robust_groupB %in% TRUE),
    `MaAsLin2` = data.frame(taxon = tx, p = mz$pval, est = mz$coef, call = mz$qval < 0.05),
    `LinDA` = data.frame(taxon = tx, p = ld$pvalue, est = ld$log2FoldChange, call = ld$padj < 0.05),
    `DESeq2` = data.frame(taxon = tx, p = ds$pvalue, est = ds$log2FoldChange, call = ds$padj < 0.05),
    `Wilcoxon CLR` = wil(clr),
    `Wilcoxon TSS` = wil(tss))
}

## per repetition and method: share of truth-0 taxa with p < 0.05, calls, any false call, power, FDP, wrong direction
rep_stats <- function(r, truth) {
  call <- !is.na(r$call) & r$call
  o <- call_outcome(call, r$est, r$taxon, truth)
  s <- call_summary(o)
  null <- truth[r$taxon] == 0
  changed <- sum(!null)
  c(p05 = mean(r$p[null] < 0.05, na.rm = TRUE), calls = s$calls, any_false = as.numeric(s$false + s$wrong_direction > 0),
    false = s$false, power = if (changed) s$right / changed else NA, fdp = s$fdp, wrong = s$wrong_direction, taxa = length(truth))
}

rows <- list()
add <- function(scenario, m, subset, metric, value, n) {
  rows[[length(rows) + 1]] <<- scores_rows(FOLDER, sprintf("smokers left throat, %d + %d samples, %d repetitions (%s)",
                                                        N_PER_GROUP, N_PER_GROUP, REPS, scenario),
                                           m, format(packageVersion(versions[[m]])), "default", subset, metric, value, n)
}
cat(sprintf("Template: smokers (left throat), %d samples with >= 1,000 reads, %d OTUs present in >= 10%%; %d OTUs eligible to change (>= 50%%); bloom OTU %s (%.1f%% of reads)\n",
            ncol(otu), nrow(otu), length(eligible), bloom_otu, 100 * sum(otu[bloom_otu, ]) / sum(otu)))
cat(sprintf("Each repetition: %d + %d samples drawn at random; OTUs in < %d samples dropped; %d repetitions per scenario, %d cores\n",
            N_PER_GROUP, N_PER_GROUP, MIN_PRESENT, REPS, CORES))
for (scenario in c("null", "planted", "bloom")) {
  res <- mclapply(seq_len(REPS), function(i) {
    d <- simulate(scenario, i)
    lapply(run_methods(d$cnt, d$grp, seed = 7000 + i), rep_stats, truth = d$truth)
  }, mc.cores = CORES)
  if (any(vapply(res, inherits, TRUE, "try-error"))) stop("a repetition failed: ", res[vapply(res, inherits, TRUE, "try-error")][[1]])
  tab <- t(sapply(METHODS, function(m) colMeans(do.call(rbind, lapply(res, `[[`, m)), na.rm = TRUE)))
  cat(sprintf("\n== %s (mean over %d repetitions; %.0f OTUs tested per repetition)\n", scenario, REPS, tab[1, "taxa"]))
  if (scenario == "null") {
    print(round(tab[, c("p05", "calls", "any_false")], 3))
    cat("   p05 = share of OTUs with p < 0.05 (expected 0.05); calls = significant OTUs per repetition (all false);\n")
    cat("   any_false = share of repetitions with at least one call (adjusted p < 0.05 keeps it <= 0.05 under the null)\n")
  } else {
    print(round(tab[, c("power", "fdp", "calls", "false", "wrong", "p05")], 3))
    cat("   power = changed OTUs called in the right direction; fdp = (false + wrong-direction calls) / calls (target 0.05);\n")
    cat("   false = unchanged OTUs called per repetition; wrong = changed OTUs called in the wrong direction; p05 = unchanged OTUs with p < 0.05\n")
  }
  for (m in METHODS) {
    v <- tab[m, ]
    add(scenario, m, "unchanged OTUs", "share_p_below_0.05", round(v[["p05"]], 4), REPS)
    add(scenario, m, "all OTUs", "calls_per_rep", round(v[["calls"]], 2), REPS)
    if (scenario == "null") {
      add(scenario, m, "all OTUs", "share_reps_any_call", round(v[["any_false"]], 4), REPS)
    } else {
      add(scenario, m, "changed OTUs", "power", round(v[["power"]], 4), REPS)
      add(scenario, m, "calls", "fdp", round(v[["fdp"]], 4), REPS)
      add(scenario, m, "unchanged OTUs", "false_calls_per_rep", round(v[["false"]], 2), REPS)
    }
  }
}
write_scores(do.call(rbind, rows), file.path(script_dir, "scores.csv"))
cat("\n", R.version.string, " | ", paste(sprintf("%s %s", unique(versions[versions != "stats"]),
    vapply(unique(versions[versions != "stats"]), function(p) format(packageVersion(p)), "")), collapse = " | "), "\n", sep = "")
