## Known-answer test of tcga_signatures.R on synthetic data in the format of the real downloads (run by CI).
##
## 500 patients, 1,500 genes. A latent proliferation factor f drives the 11 PAM50 proliferation genes (loading 0.9),
## 30% of the other genes (loading 0.3-0.8) and survival (HR 1.8 per SD of f); stage I-IV adds HR 1.5 per level and
## is independent of f. No gene affects survival except through f, so:
##   permuted labels            Cox, median split and corrected cut-point p: share p < 0.05 about 0.05;
##                              naive cut-point p: far above 0.05
##   observed labels            random sets of all genes mostly "predict" survival (they average in f-genes) ...
##                              ... but not after adjusting for the proliferation score, and sets of
##                              proliferation-free genes do not
##   proliferation signature    strongly associated
## Usage: Rscript test/test_synthetic.R   (about 30 s; writes only to a temporary directory)
suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
folder <- normalizePath(file.path(script_dir, ".."))
tmp <- tempfile("synthetic-tcga-"); dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
set.seed(7)
PROLIF <- c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C")
n <- 500; G <- 1500
genes <- c(PROLIF, sprintf("GENE%04d", seq_len(G - length(PROLIF))))
f <- rnorm(n)
load <- c(rep(0.9, length(PROLIF)), ifelse(runif(G - length(PROLIF)) < 0.3, runif(G - length(PROLIF), 0.3, 0.8), 0))
expr <- 6 + outer(load, f) + sqrt(1 - load^2) * matrix(rnorm(G * n), G)
ids <- sprintf("TCGA-SY-%04d-01", seq_len(n))
ex <- data.table(sample = genes, round(expr, 4)); setnames(ex, c("sample", ids))
fwrite(ex, file.path(tmp, "SYN_HiSeqV2.gz"), sep = "\t")
stage <- sample(1:4, n, replace = TRUE, prob = c(0.2, 0.4, 0.3, 0.1))
tt <- rexp(n, 0.1 * exp(log(1.8) * f + log(1.5) * (stage - 1))) * 365
cc <- runif(n, 0, 3650)
cdr <- data.table(sample = ids, `cancer type abbreviation` = "SYN",
                  ajcc_pathologic_tumor_stage = c("Stage I", "Stage IIA", "Stage IIIB", "Stage IV")[stage],
                  PFI = as.integer(tt <= cc), PFI.time = round(pmin(tt, cc)))
fwrite(cdr, file.path(tmp, "Survival_SupplementalTable_S1_20171025_xena_sp"), sep = "\t")

log <- system2("Rscript", c(shQuote(file.path(folder, "tcga_signatures.R")), shQuote(tmp), "SYN"), stdout = TRUE, stderr = TRUE,
               env = c("N_SETS=100", "N_PERM=100", paste0("OUT_DIR=", shQuote(tmp)), "OPENBLAS_NUM_THREADS=1"))
if (!is.null(attr(log, "status"))) { cat(log, sep = "\n"); stop("tcga_signatures.R failed") }
cat(log, sep = "\n")
p <- fread(file.path(tmp, "tcga_pvalues.csv"))
share <- function(fam, lab, a) mean(p[family == fam & labels == lab][[a]] < 0.05, na.rm = TRUE)
R <- "random 50-gene sets"; RF <- "random 50-gene sets, proliferation-free"; PR <- "proliferation (PAM50, 11 genes)"
checks <- list(
  list("permuted, Cox continuous: share p < 0.05 <= 0.12", share(R, "permuted", "Cox, continuous score") <= 0.12),
  list("permuted, median split: share p < 0.05 <= 0.12", share(R, "permuted", "median split") <= 0.12),
  list("permuted, corrected cut-point: share p < 0.05 <= 0.12", share(R, "permuted", "optimal cut-point, corrected p") <= 0.12),
  list("permuted, naive cut-point: share p < 0.05 >= 0.25", share(R, "permuted", "optimal cut-point, naive p") >= 0.25),
  list("permuted, proliferation signature, naive cut-point: share >= 0.25", share(PR, "permuted", "optimal cut-point, naive p") >= 0.25),
  list("observed, random sets, Cox continuous: share p < 0.05 >= 0.5", share(R, "observed", "Cox, continuous score") >= 0.5),
  list("observed, random sets, adjusted for proliferation: share <= 0.15", share(R, "observed", "Cox, adjusted for proliferation") <= 0.15),
  list("observed, random sets, adjusted for proliferation and stage: share <= 0.15",
       share(R, "observed", "Cox, adjusted for proliferation and stage") <= 0.15),
  list("observed, proliferation-free random sets, Cox continuous: share <= 0.15", share(RF, "observed", "Cox, continuous score") <= 0.15),
  list("observed, proliferation signature, Cox continuous: p < 1e-4", p[family == PR & labels == "observed"][["Cox, continuous score"]] < 1e-4),
  list("one row per signature and label set", nrow(p) == 2 * 2 * 100 + 2 * (1 + 100)))
cat("\n")
for (k in checks) cat(if (isTRUE(k[[2]])) "ok  " else "FAIL", k[[1]], "\n")
if (!all(vapply(checks, function(k) isTRUE(k[[2]]), logical(1)))) quit(status = 1)
cat("all synthetic checks passed\n")
