args <- commandArgs(TRUE); W <- args[1]; variant <- args[2]
if (nzchar(Sys.getenv("R_LIBS_USER"))) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(WGCNA)); enableWGCNAThreads(16); options(stringsAsFactors=FALSE)
P <- readRDS(file.path(W,"out/prep.rds")); grade <- P$grade; cpf <- as.matrix(P$cpfull)
top <- switch(variant, A=P$top, B=P$topB, C=P$topC)
X <- t(cpf[top,]); cat("== variant", variant, " dim", dim(X), " NA frac", mean(is.na(X)), "\n")
out <- list()
## 2. soft threshold
sft <- pickSoftThreshold(X, powerVector=c(1:20,22,24,26,28,30), networkType="signed", corFnc="cor", corOptions=list(use="p"), verbose=0, RsquaredCut=0.85)
fi <- sft$fitIndices; fi$signedR2 <- -sign(fi$slope)*fi$SFT.R.sq
print(round(fi[,c("Power","SFT.R.sq","slope","signedR2","mean.k.")],3))
for (cut in c(.8,.85,.9)) cat("first power signedR2 >=",cut,":", fi$Power[which(fi$signedR2>=cut)[1]], "\n")
cat("powerEstimate (cut .85):", sft$powerEstimate, "\n")
out$fi <- fi; out$pe <- sft$powerEstimate
run <- function(pw) blockwiseModules(X, power=pw, networkType="signed", TOMType="signed", minModuleSize=30, deepSplit=2, mergeCutHeight=0.25, maxBlockSize=5000, numericLabels=TRUE, pamRespectsDendro=FALSE, corType="pearson", randomSeed=12345, verbose=0)
net18 <- run(18); cat("power 18: modules", length(setdiff(unique(net18$colors),0)), " grey %", round(100*mean(net18$colors==0),1), "\n")
net <- run(10); col <- net$colors
nm <- length(setdiff(unique(col),0)); cat("power 10: modules", nm, " grey %", round(100*mean(col==0),1), "\n"); print(table(col))
ME <- net$MEs; ME <- ME[, colnames(ME)!="ME0", drop=FALSE]
r <- cor(ME, grade); p <- corPvalueStudent(r, length(grade)); rs <- cor(ME, grade, method="spearman")
mt <- data.frame(mod=colnames(ME), size=as.numeric(table(col)[sub("ME","",colnames(ME))]), r=round(r[,1],3), p=signif(p[,1],3), rho=round(rs[,1],3))
mt <- mt[order(-abs(mt$r)),]; print(mt); cat("modules p<0.05:", sum(mt$p<0.05), "\n")
topm <- as.numeric(sub("ME","",mt$mod[1]))
Xi <- impute::impute.knn(t(X))$data; Xi <- t(Xi)
kme <- cor(Xi[, col==topm], ME[,mt$mod[1]]); cat("top module hubs (kME):", paste(head(rownames(kme)[order(-kme)],15), collapse=" "), "\n")
## 4. endothelial
mk <- c("PECAM1","VWF","CDH5","ENG","EMCN","ESAM"); cat("markers in full matrix:", paste(mk[mk %in% rownames(cpf)],collapse=","), " in input set:", paste(mk[mk%in%top],collapse=","), " in top module:", paste(mk[mk %in% top[col==topm]],collapse=","), "\n")
Z <- scale(t(cpf[mk[mk %in% rownames(cpf)],])); endo <- rowMeans(Z, na.rm=TRUE)
genes_tm <- setdiff(top[col==topm], mk)
me_nomk <- moduleEigengenes(X[, genes_tm], rep(1,length(genes_tm)))$eigengenes[,1]
me_full <- ME[,mt$mod[1]]
cat("cor(ME_top, endo) full:", round(cor(me_full,endo),3), " markers removed:", round(cor(me_nomk,endo),3), "\n")
pc <- function(a,b,z) cor(resid(lm(a~z)), resid(lm(b~z)))
cat("grade r:", round(cor(me_full,grade),3), " partial|endo:", round(pc(me_full,grade,endo),3), " (markers-removed ME:", round(cor(me_nomk,grade),3), "->", round(pc(me_nomk,grade,endo),3), ")\n")
n <- length(grade); rp <- pc(me_full,grade,endo); cat("partial p ~", signif(2*pt(-abs(rp*sqrt((n-3)/(1-rp^2))), n-3),3), "\n")
cat("cor(endo, grade):", round(cor(endo,grade),3), "\n")
## 5. random sets
set.seed(1); Xs <- scale(Xi); sig <- mt[mt$p<0.05,]
pc1 <- function(M) { s <- svd(M, nu=1, nv=0); s$u[,1] }
sig$null95 <- NA; sig$exceed <- NA
for (i in seq_len(nrow(sig))) { sz <- sig$size[i]
  nr <- replicate(1000, abs(cor(pc1(Xs[, sample(ncol(Xs), sz)]), grade)))
  sig$null95[i] <- quantile(nr,.95); sig$nullmed[i] <- median(nr); sig$exceed[i] <- abs(sig$r[i]) > sig$null95[i] }
print(sig); cat("random-set: exceed 95th pct:", sum(sig$exceed), "of", nrow(sig), "\n")
out$mt <- mt; out$sig <- sig; out$col <- col; out$ME <- ME; out$endo <- endo
saveRDS(out, file.path(W, paste0("out/wgcna_",variant,".rds")))
saveRDS(list(Xi=Xi, ME=ME, grade=grade, cpf=cpf), file.path(W, paste0("out/forcv_",variant,".rds")))
