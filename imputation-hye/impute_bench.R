## Imputation methods scored against the known HYE ratios (PXD028735, MaxQuant 2.8.1 LFQ, 3 vs 3).
## Truth log2(B/A): human 0, yeast -1, E. coli +2. E. coli is 5% of sample A -> many values missing in A (MNAR).
suppressPackageStartupMessages({library(data.table); library(limma); library(imputeLCMD); library(impute)})
args <- commandArgs(TRUE); dir <- args[1]; tag <- args[2]
pg <- fread(file.path(dir, "proteinGroups.txt"), sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
sp_of <- function(h) { ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606`="Human", `OX=559292`="Yeast", `OX=83333`="Ecoli")[ox]); s <- s[!is.na(s)]
  if (length(s) == 1) s else "Other" }
pg[, species := vapply(`Fasta headers`, sp_of, "")]; pg <- pg[species %in% c("Human","Yeast","Ecoli")]
A <- paste("LFQ intensity", c("A_Alpha","A_Beta","A_Gamma")); B <- paste("LFQ intensity", c("B_Alpha","B_Beta","B_Gamma"))
m <- as.matrix(pg[, c(A, B), with = FALSE]); m[m == 0] <- NA; m <- log2(m); sp <- pg$species
truth <- c(Human = 0, Yeast = -1, Ecoli = 2)
grp <- factor(rep(c("A","B"), each = 3)); design <- model.matrix(~ grp)

perseus <- function(x, width = 0.3, shift = 1.8) { for (j in seq_len(ncol(x))) { v <- x[, j]; mu <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
  nas <- is.na(v); x[nas, j] <- rnorm(sum(nas), mu - shift * s, width * s) }; x }
methods <- list(
  none    = function(x) x,
  MinDet  = function(x) impute.MinDet(x, q = 0.01),
  MinProb = function(x) impute.MinProb(x, q = 0.01, tune.sigma = 1),
  Perseus = function(x) perseus(x),
  QRILC   = function(x) impute.QRILC(x)[[1]],
  kNN     = function(x) impute.knn(x, k = 10, rowmax = 0.9, colmax = 0.9)$data)

score <- function(filter_name, keep) {
  x0 <- m[keep, ]; s0 <- sp[keep]
  cat(sprintf("\n==== %s | filter: %s | proteins %d (human %d, yeast %d, E. coli %d) | missing cells %.1f%% ====\n",
      tag, filter_name, nrow(x0), sum(s0=="Human"), sum(s0=="Yeast"), sum(s0=="Ecoli"), 100*mean(is.na(x0))))
  cat("method   | testable | E.coli median lfc | yeast median lfc | within ±0.5 H/Y/E | changed found | human FP | FDP\n")
  out <- list()
  for (nm in names(methods)) {
    set.seed(1); x <- methods[[nm]](x0)
    ok <- rowSums(!is.na(x[, 1:3])) >= 2 & rowSums(!is.na(x[, 4:6])) >= 2
    x <- x[ok, ]; s <- s0[ok]
    lfc <- rowMeans(x[, 4:6], na.rm = TRUE) - rowMeans(x[, 1:3], na.rm = TRUE)
    lfc_c <- lfc - median(lfc[s == "Human"])          # centre on the unchanged background
    within <- sapply(c("Human","Yeast","Ecoli"), function(k) mean(abs(lfc_c[s == k] - truth[k]) <= 0.5))
    fit <- eBayes(lmFit(x, design)); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    call <- tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5
    found <- sum(call & s != "Human"); fp <- sum(call & s == "Human")
    cat(sprintf("%-8s | %8d | %17.2f | %16.2f | %3.0f/%3.0f/%3.0f %% | %5d of %5d | %8d | %4.1f%%\n", nm, nrow(x),
        median(lfc_c[s=="Ecoli"]), median(lfc_c[s=="Yeast"]), 100*within["Human"], 100*within["Yeast"], 100*within["Ecoli"],
        found, sum(s != "Human"), fp, 100 * fp / max(sum(call), 1)))
    out[[nm]] <- data.table(method = nm, filter = filter_name, data = tag, species = s, lfc_c = lfc_c, adjp = tt$adj.P.Val)
  }
  rbindlist(out)
}
both <- rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2
either <- rowSums(!is.na(m[, 1:3])) >= 2 | rowSums(!is.na(m[, 4:6])) >= 2
r1 <- score(">=2 of 3 in BOTH groups", both)
r2 <- score(">=2 of 3 in EITHER group", either)
fwrite(rbind(r1, r2), sprintf("impute_%s.csv.gz", tag))
