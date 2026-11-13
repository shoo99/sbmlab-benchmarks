## Random gene signatures and survival in TCGA (issue #11, part 2; Venet, Dumont & Detours 2011)
##
## For each cohort, every signature score (the mean z-score of its genes) is tested for association with survival:
## the four analyses of survival_tests.R, plus Cox models adjusted for proliferation, and for proliferation and stage.
## Each test runs on the observed survival labels and on permuted labels (time and event shuffled together across
## patients), where the true association is zero by construction.
##
## Signatures
##   random 50-gene sets                         N_SETS sets of 50 genes drawn from all expressed genes
##   random 50-gene sets, proliferation-free     drawn only from genes with |r| < 0.2 with the proliferation score
##   proliferation (PAM50, 11 genes)             Nielsen et al. 2010; also the adjustment covariate
##   top 50 proliferation-correlated genes       highest r with the proliferation score, the 11 genes excluded
##                                               (an in-cohort stand-in for Venet's meta-PCNA)
## The fixed signatures get N_PERM permutations each; each random set gets one permutation of its own.
##
## Data (download.sh): Xena TCGA hub HiSeqV2 (RSEM, log2(norm + 1)) per cohort, and the curated endpoints of the
## TCGA Pan-Cancer Clinical Data Resource (Liu et al. 2018; Xena: Survival_SupplementalTable_S1_20171025_xena_sp).
## Patients: primary tumours (sample code 01) with the endpoint available and time > 0. Genes: above 0 in at least
## half of the tumours and not constant; z-scored within the cohort. Stage: AJCC pathologic stage I-IV (sub-stages
## merged; stage X and missing values dropped from the stage-adjusted model only).
##
## Usage: Rscript tcga_signatures.R [data_dir] [cohort ...] > tcga_signatures.out
##   defaults: data_dir = data next to this script, cohorts BRCA KIRC LUAD COAD
##   env: ENDPOINT (PFI), N_SETS (1000), N_PERM (200), OUT_DIR (this folder)
## Writes tcga_pvalues.csv (one row per signature, label set and cohort; one column per analysis); score.R turns it
## into scores.csv. Seed: set.seed(11) before each cohort.
suppressPackageStartupMessages({ library(survival); library(maxstat); library(data.table) })
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "survival_tests.R"))
args <- commandArgs(TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path(script_dir, "data")
cohorts <- if (length(args) >= 2) args[-1] else c("BRCA", "KIRC", "LUAD", "COAD")
ENDPOINT <- Sys.getenv("ENDPOINT", "PFI")
N_SETS <- as.integer(Sys.getenv("N_SETS", "1000"))
N_PERM <- as.integer(Sys.getenv("N_PERM", "200"))
OUT_DIR <- Sys.getenv("OUT_DIR", script_dir)
SET_SIZE <- 50
PROLIF <- c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C")
ALL_ANALYSES <- c(ANALYSES, "Cox, adjusted for proliferation", "Cox, adjusted for proliferation and stage")
SHORT <- c("cox", "median", "opt_naive", "opt_corrected", "cox_adj_prolif", "cox_adj_prolif_stage")

surv_file <- file.path(data_dir, "Survival_SupplementalTable_S1_20171025_xena_sp")
cdr <- fread(surv_file, sep = "\t", colClasses = "character", na.strings = c("", "NA"))
stage_of <- function(s) {
  r <- regmatches(s, regexpr("^Stage (IV|I{1,3})", s))
  out <- rep(NA_character_, length(s)); out[grepl("^Stage (IV|I{1,3})", s)] <- sub("Stage ", "", r)
  out
}

run_tests <- function(time, event, score, prolif, stage) {
  a <- marker_tests(time, event, score)$p
  ok <- !is.na(stage)
  adj1 <- if (is.null(prolif)) NA else cox_p(time, event, score, data.frame(prolif = prolif))
  adj2 <- if (is.null(prolif) || length(unique(stage[ok])) < 2) NA else
    cox_p(time[ok], event[ok], score[ok], data.frame(prolif = prolif[ok], stage = factor(stage[ok])))
  setNames(c(a, adj1, adj2), ALL_ANALYSES)
}

out <- list()
for (cohort in cohorts) {
  set.seed(11)
  f <- file.path(data_dir, paste0(cohort, "_HiSeqV2.gz"))
  ex <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", data.table = FALSE)
  genes <- ex[[1]]
  m <- as.matrix(ex[, -1]); rownames(m) <- genes
  m <- m[, grepl("^TCGA-..-....-01$", colnames(m)), drop = FALSE]
  cl <- cdr[cdr$sample %in% colnames(m), ]
  time <- as.numeric(cl[[paste0(ENDPOINT, ".time")]]); event <- as.integer(cl[[ENDPOINT]])
  keep <- !is.na(time) & !is.na(event) & time > 0
  cl <- cl[keep, ]; time <- time[keep]; event <- event[keep]
  m <- m[, cl$sample, drop = FALSE]
  stage <- stage_of(cl$ajcc_pathologic_tumor_stage)
  expressed <- rowMeans(m > 0) >= 0.5 & apply(m, 1, sd) > 0 & !duplicated(rownames(m))
  m <- m[expressed, , drop = FALSE]
  z <- (m - rowMeans(m)) / apply(m, 1, sd)
  score_of <- function(g) colMeans(z[g, , drop = FALSE])

  prolif_genes <- intersect(PROLIF, rownames(z))
  prolif <- score_of(prolif_genes)
  r <- as.vector(cor(t(z), prolif))
  names(r) <- rownames(z)
  top_prolif <- names(sort(r[!names(r) %in% PROLIF], decreasing = TRUE))[seq_len(SET_SIZE)]
  free_pool <- names(r)[abs(r) < 0.2]
  cat(sprintf("== %s: %d patients, %d %s events (%.0f%%), stage known for %d (%s); %d expressed genes\n",
              cohort, length(time), sum(event), ENDPOINT, 100 * mean(event), sum(!is.na(stage)),
              paste(names(table(stage)), table(stage), sep = ":", collapse = " "), nrow(z)))
  cat(sprintf("   proliferation score from %d of %d genes; %d genes with |r| < 0.2 (proliferation-free pool), median |r| of all genes %.2f\n",
              length(prolif_genes), length(PROLIF), length(free_pool), median(abs(r))))
  base <- data.frame(cohort = cohort, endpoint = ENDPOINT, patients = length(time), events = sum(event))
  record <- function(family, signature, labels, p) {
    out[[length(out) + 1]] <<- cbind(base, family = family, signature = signature, labels = labels,
                                    as.data.frame(t(signif(p, 4)), check.names = FALSE))
  }
  # random families: observed labels, then one permutation per set
  for (fam in list(list("random 50-gene sets", rownames(z)), list("random 50-gene sets, proliferation-free", free_pool))) {
    for (s in seq_len(N_SETS)) {
      sc <- score_of(sample(fam[[2]], SET_SIZE))
      record(fam[[1]], paste0("set", s), "observed", run_tests(time, event, sc, prolif, stage))
      pi <- sample(length(time))
      record(fam[[1]], paste0("set", s), "permuted", run_tests(time[pi], event[pi], sc, prolif, stage))
    }
  }
  # fixed signatures: observed labels and N_PERM permutations (proliferation is not adjusted for itself)
  for (fx in list(list("proliferation (PAM50, 11 genes)", prolif_genes, FALSE),
                  list("top 50 proliferation-correlated genes", top_prolif, TRUE))) {
    sc <- score_of(fx[[2]]); pr <- if (fx[[3]]) prolif else NULL
    record(fx[[1]], fx[[1]], "observed", run_tests(time, event, sc, pr, stage))
    for (k in seq_len(N_PERM)) {
      pi <- sample(length(time))
      record(fx[[1]], paste0("perm", k), "permuted", run_tests(time[pi], event[pi], sc, pr, stage))
    }
  }
  res <- do.call(rbind, out[vapply(out, function(x) x$cohort[1] == cohort, logical(1))])
  obs <- res[res$labels == "observed" & res$family == "random 50-gene sets", ALL_ANALYSES]
  cat("   random 50-gene sets, observed labels, share p < 0.05:\n")
  print(setNames(round(colMeans(obs < 0.05, na.rm = TRUE), 3), SHORT))
}
res <- do.call(rbind, out)
fwrite(res, file.path(OUT_DIR, "tcga_pvalues.csv"))
cat(sprintf("\nwrote tcga_pvalues.csv: %d rows (%s)\n", nrow(res), paste(cohorts, collapse = ", ")))
cat(R.version.string, " | survival ", format(packageVersion("survival")), " | maxstat ", format(packageVersion("maxstat")),
    " | data.table ", format(packageVersion("data.table")), "\n", sep = "")
