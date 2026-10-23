args <- commandArgs(TRUE); W <- args[1]
if (nzchar(Sys.getenv("R_LIBS_USER"))) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
cp <- read.delim(file.path(W,"data/HS_CPTAC_CCRCC_proteome_Tumor.cct"), row.names=1, check.names=FALSE)
cli <- read.delim(file.path(W,"data/HS_CPTAC_CCRCC_CLI.tsi"), check.names=FALSE)[-1,]
cat("CPTAC proteome:", dim(cp), " NA frac total:", mean(is.na(cp)), "\n")
cli <- cli[cli$Case_ID %in% colnames(cp),]
cat("clinical rows matched:", nrow(cli), "\n"); print(table(cli$Histologic_Grade, useNA="a"))
cat("all clinical rows grade table:\n"); cli0 <- read.delim(file.path(W,"data/HS_CPTAC_CCRCC_CLI.tsi"))[-1,]; print(table(cli0$Histologic_Grade,useNA="a"))
tj <- read.csv(file.path(W,"data/tjrcc/matrix_protein.csv"), row.names=1, check.names=FALSE)
tj <- tj[, grepl("_T$", colnames(tj))]; tj[tj==0] <- NA
cat("TJ tumours:", ncol(tj), " proteins:", nrow(tj), " NA frac:", mean(is.na(tj)), " zeros:", mean(tj==0,na.rm=TRUE), "\n")
print(summary(as.numeric(as.matrix(tj))))
mc <- rowMeans(is.na(cp)); mt <- rowMeans(is.na(tj))
shared <- intersect(rownames(cp), rownames(tj))
keep <- shared[mc[shared] <= 0.10 & mt[shared] <= 0.10]
keepC <- rownames(cp)[mc <= 0.10]
cat("shared:", length(shared), " shared & <=10% missing both:", length(keep), " CPTAC-only <=10%:", length(keepC), "\n")
v <- apply(cp[keep,], 1, var, na.rm=TRUE)
top <- names(sort(v, decreasing=TRUE))[1:3000]
vC <- apply(cp[keepC,],1,var,na.rm=TRUE); topC <- names(sort(vC,decreasing=TRUE))[1:3000]
cat("overlap top3000 shared vs CPTAC-only:", length(intersect(top,topC)), "\n")
X <- t(as.matrix(cp[top, cli$Case_ID]))
cat("NA frac in X:", mean(is.na(X)), "\n")
grade <- as.numeric(sub("G","",cli$Histologic_Grade))
tj0 <- read.csv(file.path(W,"data/tjrcc/matrix_protein.csv"), row.names=1, check.names=FALSE); tj0 <- tj0[,grepl("_T$",colnames(tj0))]
keepB <- shared[mc[shared] <= 0.10 & rowMeans(is.na(tj0[shared,]))<=0.10]; vB <- apply(cp[keepB,],1,var,na.rm=TRUE); topB <- names(sort(vB,decreasing=TRUE))[1:3000]
cat("variant B (zeros as values) n:", length(keepB), " overlap with A:", length(intersect(top,topB)), "\n")
saveRDS(list(topB=topB, X=X, grade=grade, ids=cli$Case_ID, cpfull=cp[, cli$Case_ID], tj=tj, top=top, topC=topC), file.path(W,"out/prep.rds"))
