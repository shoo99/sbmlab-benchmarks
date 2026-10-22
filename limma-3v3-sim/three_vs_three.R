# Numbers behind "Three samples per group: switching from the t-test to limma" (sbmlab.com)
# Rscript three_vs_three.R   (R 4.x, limma 3.58.1; ~3 min)
suppressMessages(library(limma))

cat("== 1. Smallest attainable two-sided Mann-Whitney p (no ties)\n")
for (n in 3:5) cat(sprintf("n=%d per group: %d splits, min p = %.4f\n", n, choose(2*n, n), 2/choose(2*n, n)))

cat("\n== 2. Power of a two-sided two-sample t-test, alpha 0.05, n = 3 per group\n")
for (d in c(1, 2, 3)) cat(sprintf("d = %g: power %.3f\n", d, power.t.test(n=3, delta=d, sd=1, sig.level=0.05)$power))
for (d in c(1, 1.5)) cat(sprintf("d = %g: n per group for 80%% power = %.2f\n", d, power.t.test(power=0.8, delta=d, sd=1, sig.level=0.05)$n))

cat("\n== 3. Synthetic experiment with a known answer\n")
# 400 features, 40 truly changed (log2 +1.5 / -1.5 alternating), 3 vs 3,
# per-feature within-group SD ~ log-normal(median 0.5 on log2 scale, sdlog 0.6), 500 repetitions.
# t-test = Student (equal variances); limma = lmFit + eBayes defaults; fold change = difference of group means (log2).
set.seed(1); R <- 500; P <- 400; T <- 40; n <- 3
grp <- factor(rep(c("A", "B"), each=n)); design <- model.matrix(~grp)
rules <- c("p<0.05, t-test", "p<0.05, limma", "p<0.05 & >=2-fold, t-test", "p<0.05 & >=2-fold, limma", "BH 5%, t-test", "BH 5%, limma")
acc <- matrix(0, 6, 3, dimnames=list(rules, c("selected", "false", "true")))
for (r in 1:R) {
  sd <- rlnorm(P, log(0.5), 0.6)
  eff <- c(rep(c(1.5, -1.5), length.out=T), rep(0, P - T)); truth <- eff != 0
  X <- matrix(rnorm(P * 2 * n, 0, rep(sd, 2 * n)), P); X[, (n+1):(2*n)] <- X[, (n+1):(2*n)] + eff
  lfc <- rowMeans(X[, (n+1):(2*n)]) - rowMeans(X[, 1:n])
  pt <- apply(X, 1, function(x) t.test(x[(n+1):(2*n)], x[1:n], var.equal=TRUE)$p.value)
  pl <- eBayes(lmFit(X, design))$p.value[, 2]
  sets <- list(pt < 0.05, pl < 0.05, pt < 0.05 & abs(lfc) >= 1, pl < 0.05 & abs(lfc) >= 1,
               p.adjust(pt, "BH") < 0.05, p.adjust(pl, "BH") < 0.05)
  for (k in 1:6) acc[k, ] <- acc[k, ] + c(sum(sets[[k]]), sum(sets[[k]] & !truth), sum(sets[[k]] & truth))
}
acc <- acc / R
print(cbind(round(acc, 1), false_share_pct=round(100 * acc[, 2] / acc[, 1])))
cat(sprintf("\nExpected false hits at p<0.05 if the test were exact: %d nulls x 0.05 = %.1f\n", P - T, (P - T) * 0.05))
sessionInfo()$R.version$version.string
