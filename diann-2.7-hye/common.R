suppressPackageStartupMessages(library(data.table))
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)
sp_of <- function(n){s<-unique(regmatches(n,gregexpr("_(HUMAN|YEAST|ECOLI)\\b",n))[[1]]);s<-c(`_HUMAN`="Human",`_YEAST`="Yeast",`_ECOLI`="Ecoli")[s];if(length(s)==1)s else NA_character_}
fcs <- function(f){x<-fread(f,sep="\t");rc<-grep("Condition_",names(x),value=TRUE);A<-grep("Condition_A",rc,value=TRUE);B<-grep("Condition_B",rc,value=TRUE)
 m<-as.matrix(x[,c(A,B),with=FALSE]);m[m==0]<-NA;m<-log2(m);ok<-rowSums(!is.na(m[,A]))>=2&rowSums(!is.na(m[,B]))>=2
 data.table(pg=x$Protein.Group,sp=vapply(x$Protein.Names,sp_of,""),fc=rowMeans(m[,B],na.rm=TRUE)-rowMeans(m[,A],na.rm=TRUE),ok=ok,
            intA=rowMeans(m[,A],na.rm=TRUE))}
a<-fcs("../diann-hye/report.pg_matrix.tsv")[ok==TRUE]; b<-fcs("report.pg_matrix.tsv")[ok==TRUE]
j<-merge(a,b,by=c("pg","sp"),suffixes=c(".232",".270"))
cat("common quantified:",nrow(j)," only 2.3.2:",nrow(a)-nrow(j)," only 2.7.0:",nrow(b)-nrow(j),"\n")
for(s in names(EXPECT)){e<-EXPECT[[s]];q<-j[sp==s]
 cat(sprintf("%s common n=%d  within0.5 2.3.2=%.3f 2.7.0=%.3f  median %.3f / %.3f  IQR %.3f / %.3f\n",s,nrow(q),mean(abs(q$fc.232-e)<=.5),mean(abs(q$fc.270-e)<=.5),median(q$fc.232),median(q$fc.270),IQR(q$fc.232),IQR(q$fc.270)))
 o<-b[sp==s & !(pg %in% a$pg)]
 if(nrow(o)) cat(sprintf("   %s only-2.7.0 n=%d within0.5=%.3f median %.3f  median intA %.2f (common %.2f)\n",s,nrow(o),mean(abs(o$fc-e)<=.5),median(o$fc),median(o$intA,na.rm=TRUE),median(q$intA.270,na.rm=TRUE)))}
## Of the protein groups only 2.7.0 quantified, how many were already in the 2.3.2 matrix (with too few values to pass 2 of 3)?
x <- fread("../diann-hye/report.pg_matrix.tsv"); nw <- b[!(pg %in% a$pg)]
cat("\nonly-2.7.0 groups already present in the 2.3.2 matrix:\n"); print(table(species = nw$sp, in_2.3.2_matrix = nw$pg %in% x$Protein.Group))
