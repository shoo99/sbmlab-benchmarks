lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, showWarnings=FALSE, recursive=TRUE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=16)
need <- c("impute","preprocessCore")
miss <- need[!sapply(need, requireNamespace, quietly=TRUE)]
if (length(miss)) { if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", lib=lib); BiocManager::install(miss, lib=lib, update=FALSE, ask=FALSE) }
install.packages("WGCNA", lib=lib, dependencies=c("Depends","Imports"))
library(WGCNA); cat("WGCNA", as.character(packageVersion("WGCNA")), "\n")
