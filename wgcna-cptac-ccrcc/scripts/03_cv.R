args <- commandArgs(TRUE); W <- args[1]; v <- args[2]
F <- readRDS(file.path(W, paste0("out/forcv_",v,".rds"))); Xi <- F$Xi; ME <- as.matrix(F$ME); g <- F$grade
y <- as.integer(g >= 3); cat("variant", v, " G3-4:", sum(y), " G1-2:", sum(!y), " nME:", ncol(ME), "\n")
auc <- function(s, y) { r <- rank(s); n1 <- sum(y); n0 <- sum(!y); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
fitpred <- function(a, b, yy) { d <- data.frame(y=yy, a); suppressWarnings(m <- glm(y~., data=d, family=binomial)); suppressWarnings(predict(m, newdata=data.frame(b), type="link")) }
score50 <- function(Xtr, gtr, Xte, k=50) { r <- cor(Xtr, gtr)[,1]; sel <- order(-abs(r))[1:k]; mu <- colMeans(Xtr[,sel]); sd <- apply(Xtr[,sel],2,sd)
  f <- function(M) as.matrix(sweep(sweep(M[,sel,drop=FALSE],2,mu),2,sd,"/")) %*% sign(r[sel]) / k; list(tr=f(Xtr), te=f(Xte)) }
selFull <- order(-abs(cor(Xi, g)[,1]))[1:50]
R <- 50; res <- matrix(NA, R, 6, dimnames=list(NULL, c("top50_inside","top50_fulldata_sel","top50_inside_binarysel","PC15","ME_all","ME15_or_all")))
set.seed(2026)
for (rep in 1:R) { fold <- integer(length(y)); for (cl in 0:1) { idx <- which(y==cl); fold[idx] <- sample(rep_len(1:5, length(idx))) }
  P <- matrix(NA, length(y), 6)
  for (k in 1:5) { tr <- fold!=k; te <- fold==k
    s <- score50(Xi[tr,], g[tr], Xi[te,]); P[te,1] <- fitpred(data.frame(s=s$tr), data.frame(s=s$te), y[tr])
    mu <- colMeans(Xi[tr,selFull]); sdv <- apply(Xi[tr,selFull],2,sd); rr <- sign(cor(Xi[,selFull], g)[,1])
    sf <- function(M) as.matrix(sweep(sweep(M[,selFull],2,mu),2,sdv,"/")) %*% rr / 50
    P[te,2] <- fitpred(data.frame(s=sf(Xi[tr,])), data.frame(s=sf(Xi[te,])), y[tr])
    s3 <- score50(Xi[tr,], y[tr], Xi[te,]); P[te,3] <- fitpred(data.frame(s=s3$tr), data.frame(s=s3$te), y[tr])
    pcs <- prcomp(Xi[tr,], center=TRUE, scale.=TRUE); Ptr <- pcs$x[,1:15]; Pte <- predict(pcs, Xi[te,])[,1:15]
    P[te,4] <- fitpred(data.frame(Ptr), data.frame(Pte), y[tr])
    P[te,5] <- fitpred(data.frame(ME[tr,]), data.frame(ME[te,]), y[tr])
  }
  res[rep,1:5] <- apply(P[,1:5], 2, auc, y=y) }
res <- res[,1:5]
print(round(rbind(mean=colMeans(res), sd=apply(res,2,sd)),3))
cat("in-sample (no CV) AUC of full-data top50 score:", round(auc(as.matrix(scale(Xi[,selFull])) %*% sign(cor(Xi[,selFull],g)[,1]), y),3), "\n")
