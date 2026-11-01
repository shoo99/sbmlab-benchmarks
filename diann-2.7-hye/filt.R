suppressPackageStartupMessages(library(data.table))
EXPECT <- c(Human = 0, Yeast = -1, Ecoli = 2)
sp_of <- function(n){s<-unique(regmatches(n,gregexpr("_(HUMAN|YEAST|ECOLI)\\b",n))[[1]]);s<-c(`_HUMAN`="Human",`_YEAST`="Yeast",`_ECOLI`="Ecoli")[s];if(length(s)==1)s else NA_character_}
run <- function(f,lab){x<-fread(f,sep="\t");rc<-grep("Condition_",names(x),value=TRUE);A<-grep("Condition_A",rc,value=TRUE);B<-grep("Condition_B",rc,value=TRUE)
 m<-as.matrix(x[,c(A,B),with=FALSE]);m[m==0]<-NA;m<-log2(m);sp<-vapply(x$Protein.Names,sp_of,"")
 fc<-rowMeans(m[,B],na.rm=TRUE)-rowMeans(m[,A],na.rm=TRUE)
 for(rule in c("2of3","3of3")){k<-if(rule=="2of3")2 else 3
  ok<-rowSums(!is.na(m[,A]))>=k&rowSums(!is.na(m[,B]))>=k
  cat(lab,rule,paste(sapply(names(EXPECT),function(s){i<-ok&sp%in%s;sprintf("%s n=%d %.1f%%",s,sum(i),100*mean(abs(fc[i]-EXPECT[[s]])<=.5))}),collapse=" | "),"\n")}}
run("../diann-hye/report.pg_matrix.tsv","2.3.2")
run("report.pg_matrix.tsv","2.7.0")
