suppressPackageStartupMessages(library(data.table))
A <- c("A_Alpha","A_Beta","A_Gamma"); B <- c("B_Alpha","B_Beta","B_Gamma"); EXPECT <- c(Human=0, Yeast=-1, Ecoli=2)
fcs <- function(m) { m[m==0] <- NA; m <- log2(m); colnames(m) <- c(A,B)
  ok <- rowSums(!is.na(m[,A]))>=2 & rowSums(!is.na(m[,B]))>=2
  list(ok=ok, fc=rowMeans(m[,B],na.rm=TRUE)-rowMeans(m[,A],na.rm=TRUE)) }
mq <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep="\t", quote="")
mq <- mq[!(Decoy %in% "+") & !(`Potential contaminant` %in% "+") & !(`Only identified by site` %in% "+")]
mq[, acc := sub(";.*","",`Majority protein IDs`)]
r <- fcs(as.matrix(mq[, paste("LFQ intensity", c(A,B)), with=FALSE])); mq[, `:=`(ok_mq=r$ok, fc_mq=r$fc)]
fp <- fread("work/combined_protein.tsv", sep="\t", quote="")
fp[, acc := `Protein ID`]
r <- fcs(as.matrix(fp[, paste(c(A,B), "MaxLFQ Intensity"), with=FALSE])); fp[, `:=`(ok_fp=r$ok, fc_fp=r$fc)]
fp[, species := fifelse(grepl("_HUMAN", Protein),"Human", fifelse(grepl("_YEAST",Protein),"Yeast", fifelse(grepl("_ECOLI",Protein),"Ecoli",NA_character_)))]
both <- merge(mq[ok_mq==TRUE, .(acc, fc_mq)], fp[ok_fp==TRUE & !is.na(species), .(acc, species, fc_fp)], by="acc")
cat("proteins scoreable in both (matched on leading accession):", nrow(both), "\n")
print(both[, .(n=.N, within_mq=round(mean(abs(fc_mq-EXPECT[species])<=0.5),3), within_fp=round(mean(abs(fc_fp-EXPECT[species])<=0.5),3),
               IQR_mq=round(IQR(fc_mq),3), IQR_fp=round(IQR(fc_fp),3), med_mq=round(median(fc_mq),3), med_fp=round(median(fc_fp),3)), by=species])
only_fp <- fp[ok_fp==TRUE & !is.na(species) & !(acc %in% mq$acc[mq$ok_mq])]
cat("scoreable only in FragPipe:", nrow(only_fp), "\n")
print(only_fp[, .(n=.N, within=round(mean(abs(fc_fp-EXPECT[species])<=0.5),3)), by=species])
## centred on each tool's own human median (removes the normalisation offset)
hm <- both[species=="Human", .(m_mq=median(fc_mq), m_fp=median(fc_fp))]
both[, `:=`(c_mq=fc_mq-hm$m_mq, c_fp=fc_fp-hm$m_fp)]
cat("\n== centred on human median ==\n")
print(both[, .(n=.N, within_mq=round(mean(abs(c_mq-EXPECT[species])<=0.5),3), within_fp=round(mean(abs(c_fp-EXPECT[species])<=0.5),3),
               med_mq=round(median(c_mq),3), med_fp=round(median(c_fp),3)), by=species])
