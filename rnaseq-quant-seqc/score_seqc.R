## Score each quantifier's gene matrix against the SEQC known answers, with the shared rules in common/known_ratio.R.
##   1. ERCC: log2(A/B) per spike vs log2(Mix 1 / Mix 2) of its subgroup (A carries Mix 1, B Mix 2); share within
##      +-0.5, uncentred and centred on the 1:1 subgroup (which absorbs the mRNA-content difference between A and B).
##   2. Titration order: endogenous genes with |log2 A/B| >= 1 should satisfy A > C > D > B (or the reverse).
##   3. Mixture model: C and D are mixtures of A and B by total RNA (3:1, 1:3); fit w = mRNA content of A relative to B
##      and report the median absolute log2 residual of C and D.
## Expression = counts per million endogenous reads; zeros are missing; a sample's mean needs >= min_valid replicates.
## Usage: Rscript score_seqc.R <dir with genes_*.tsv.gz> <ercc_truth.tsv> [out dir = .] [dataset label] [min_valid = 3]
## Prints the results (run.sh in the README saves them as score_seqc.out) and writes scores.csv to <out dir>.
suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "common", "known_ratio.R"))
a <- commandArgs(TRUE)
IN <- a[1]; TRUTH <- a[2]
OUT <- if (length(a) > 2) a[3] else "."
DATASET <- if (length(a) > 3) a[4] else "SEQC"
MIN_VALID <- if (length(a) > 4) as.integer(a[5]) else 3L
FOLDER <- "rnaseq-quant-seqc"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ercc <- fread(TRUTH)
ercc[, label := sprintf("ERCC %s (Mix1/Mix2 = %s)", subgroup, format(round(2^log2_mix1_over_mix2, 2)))]
truth <- setNames(ercc[, .(v = log2_mix1_over_mix2[1]), by = label]$v, ercc[, .(v = log2_mix1_over_mix2[1]), by = label]$label)
centre_label <- names(truth)[abs(truth) < 1e-9]
versions <- if (file.exists(file.path(IN, "versions.tsv"))) fread(file.path(IN, "versions.tsv")) else data.table(tool = character(), version = character())
tool_names <- c(salmon = "salmon", kallisto = "kallisto", rsem = "RSEM", featurecounts = "featureCounts")
f_from_A <- function(share_A, w) share_A * w / (share_A * w + (1 - share_A))     # mRNA fraction from A

rows <- list()
for (tool in names(tool_names)) {
  f <- file.path(IN, sprintf("genes_%s.tsv.gz", tool))
  if (!file.exists(f)) next
  x <- fread(f)
  smp <- setdiff(names(x), "gene_id")
  cols <- lapply(c(A = "A", B = "B", C = "C", D = "D"), function(s) grep(sprintf("^%s_", s), smp, value = TRUE))
  m <- as.matrix(x[, ..smp])
  is_ercc <- grepl("^ERCC-", x$gene_id)
  cpm <- sweep(m, 2, colSums(m[!is_ercc, , drop = FALSE]), "/") * 1e6
  version <- versions$version[match(tool_names[[tool]], versions$tool)]
  if (!length(version) || is.na(version)) version <- ""
  cat(sprintf("\n==== %s %s | %d genes (%d ERCC) | samples A/B/C/D: %s\n", tool_names[[tool]], version, nrow(x), sum(is_ercc),
              paste(lengths(cols), collapse = "/")))

  ## 1. ERCC, log2(A/B)
  e <- cpm[is_ercc, , drop = FALSE]
  lab <- ercc$label[match(x$gene_id[is_ercc], ercc$ercc_id)]
  r <- log2_ratio(e, a = cols$B, b = cols$A, min_valid = MIN_VALID)
  s <- score_groups(r$fc, r$ok, lab, truth, centre = centre_label)
  print(s, row.names = FALSE)
  rows[[length(rows) + 1]] <- rbind(
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", s$group, "median_log2_A_over_B", s$median, s$n),
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", s$group, "within_0.5", s$within, s$n),
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", s$group, "within_0.5_centred", s$within_centred, s$n))

  ## 2. titration order on endogenous genes
  g <- group_log2_means(cpm[!is_ercc, , drop = FALSE], cols, min_valid = MIN_VALID)
  g <- g[stats::complete.cases(g), , drop = FALSE]
  inf <- abs(g[, "A"] - g[, "B"]) >= 1
  gi <- g[inf, , drop = FALSE]
  consistent <- (gi[, "A"] > gi[, "C"] & gi[, "C"] > gi[, "D"] & gi[, "D"] > gi[, "B"]) |
                (gi[, "A"] < gi[, "C"] & gi[, "C"] < gi[, "D"] & gi[, "D"] < gi[, "B"])
  share <- round(mean(consistent), 3)
  cat(sprintf("titration order A>C>D>B (or reverse): %.3f of %d genes with |log2 A/B| >= 1 (%d genes with all four samples)\n",
              share, nrow(gi), nrow(g)))

  ## 3. mixture model: fit w on the informative genes, residuals of C and D
  la <- 2^gi[, "A"]; lb <- 2^gi[, "B"]
  resid <- function(w) {
    fc <- f_from_A(0.75, w); fd <- f_from_A(0.25, w)
    c(gi[, "C"] - log2(fc * la + (1 - fc) * lb), gi[, "D"] - log2(fd * la + (1 - fd) * lb))
  }
  w_hat <- exp(optimize(function(lw) sum(resid(exp(lw))^2), c(log(0.1), log(10)))$minimum)
  rr <- resid(w_hat); n_inf <- nrow(gi)
  mad_c <- round(median(abs(rr[seq_len(n_inf)])), 3); mad_d <- round(median(abs(rr[n_inf + seq_len(n_inf)])), 3)
  cat(sprintf("mixture model: w (mRNA content of A relative to B) = %.3f; median |log2 residual| C %.3f, D %.3f\n",
              w_hat, mad_c, mad_d))
  sub_end <- "endogenous, |log2 A/B| >= 1"
  rows[[length(rows) + 1]] <- rbind(
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", sub_end, "titration_order_consistent", share, nrow(gi)),
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", sub_end, "mixture_w", round(w_hat, 3), nrow(gi)),
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", sub_end, "mixture_median_abs_resid_C", mad_c, nrow(gi)),
    scores_rows(FOLDER, DATASET, tool_names[[tool]], version, "default", sub_end, "mixture_median_abs_resid_D", mad_d, nrow(gi)))
}
if (!length(rows)) stop("no genes_<tool>.tsv.gz in ", IN)
write_scores(do.call(rbind, rows), file.path(OUT, "scores.csv"))
cat(sprintf("\nmin_valid = %d replicates per sample; scores.csv written to %s\n", MIN_VALID, OUT))
cat(R.version.string, "| data.table", format(packageVersion("data.table")), "\n")
