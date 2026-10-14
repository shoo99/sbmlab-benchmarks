suppressPackageStartupMessages(library(data.table))
pg <- fread("/home/sysoft/mq-hye/raw/combined/txt/proteinGroups.txt", sep="\t", quote="")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
m <- as.matrix(pg[, paste("LFQ intensity", c("A_Alpha","A_Beta","A_Gamma","B_Alpha","B_Beta","B_Gamma")), with=FALSE]); m[m==0] <- NA; m <- log2(m)
library(limma); library(imputeLCMD); library(impute)

# m: proteins x samples, log2 LFQ with NA for missing; groups: A A A B B B
keep_strict <- rowSums(!is.na(m[, 1:3])) >= 2 & rowSums(!is.na(m[, 4:6])) >= 2
keep_loose  <- rowSums(!is.na(m[, 1:3])) >= 2 | rowSums(!is.na(m[, 4:6])) >= 2

perseus <- function(x, width = 0.3, shift = 1.8) {    # per sample, like Perseus' defaults
  for (j in seq_len(ncol(x))) {
    v <- x[, j]; mu <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
    x[is.na(v), j] <- rnorm(sum(is.na(v)), mu - shift * s, width * s)
  }
  x
}

set.seed(1)
x <- m[keep_loose, ]
x_mindet  <- impute.MinDet(x, q = 0.01)
x_minprob <- impute.MinProb(x, q = 0.01, tune.sigma = 1)
x_perseus <- perseus(x)
x_qrilc   <- impute.QRILC(x)[[1]]
x_knn     <- impute.knn(x, k = 10, rowmax = 0.9, colmax = 0.9)$data   # see finding 3 before using

grp <- factor(rep(c("A", "B"), each = 3))
fit <- eBayes(lmFit(x_mindet, model.matrix(~ grp)))
res <- topTable(fit, coef = 2, number = Inf)

cat("rows", nrow(x), "sig", sum(res$adj.P.Val < 0.05), "\n")
