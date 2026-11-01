## Per-protein log2(B/A) and mean log2 intensity in condition A for 2.3.2 and 2.7.0 (rule: >=2 of 3 per condition), for the figure.
suppressPackageStartupMessages(library(data.table))
sp_of <- function(n){s<-unique(regmatches(n,gregexpr("_(HUMAN|YEAST|ECOLI)\\b",n))[[1]]);s<-c(`_HUMAN`="Human",`_YEAST`="Yeast",`_ECOLI`="Ecoli")[s];if(length(s)==1)s else NA_character_}
fcs <- function(f){x<-fread(f,sep="\t");rc<-grep("Condition_",names(x),value=TRUE);A<-grep("Condition_A",rc,value=TRUE);B<-grep("Condition_B",rc,value=TRUE)
 m<-as.matrix(x[,c(A,B),with=FALSE]);m[m==0]<-NA;m<-log2(m);ok<-rowSums(!is.na(m[,A]))>=2&rowSums(!is.na(m[,B]))>=2
 data.table(pg=x$Protein.Group,species=vapply(x$Protein.Names,sp_of,""),fc=rowMeans(m[,B],na.rm=TRUE)-rowMeans(m[,A],na.rm=TRUE),intA=rowMeans(m[,A],na.rm=TRUE))[ok]}
a<-fcs(file.path("..","diann-hye","report.pg_matrix.tsv")); b<-fcs("report.pg_matrix.tsv")
b[, set := ifelse(pg %in% a$pg, "both", "only_2.7.0")]; a[, set := ifelse(pg %in% b$pg, "both", "only_2.3.2")]
fwrite(rbind(cbind(version="2.3.2",a),cbind(version="2.7.0",b)),"fc_both_versions.csv")
