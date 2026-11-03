## Read MaxQuant evidence.txt (146 MB) + aggregate per protein; median of 5 runs.
suppressPackageStartupMessages(library(data.table))
f <- commandArgs(TRUE)[1]
cols <- c("Leading razor proteins", "PEP", "Intensity", "Decoy")
t_dt <- replicate(5, system.time({
  d <- fread(f, sep = "\t", quote = "", select = cols)
  a <- d[!(Decoy %in% "+") & PEP < 0.01, .(mean_int = mean(Intensity, na.rm = TRUE), n = .N), by = `Leading razor proteins`]
})["elapsed"])
## base R as written in the post: read.delim parses and type-converts every column (several dozen), while fread,
## pandas and polars read only the 4 columns above; aggregate() also drops NA-Intensity rows (na.action = na.omit).
t_base <- replicate(5, system.time({
  d <- read.delim(f, quote = "", check.names = FALSE)
  d <- d[!(d$Decoy %in% "+") & d$PEP < 0.01, ]
  a <- aggregate(Intensity ~ `Leading razor proteins`, data = d, FUN = function(x) mean(x, na.rm = TRUE))
})["elapsed"])
## base R doing the same work as the other three (added 2026-11-02): only the 4 columns, NA rows kept, mean and n
hdr <- names(read.delim(f, quote = "", check.names = FALSE, nrows = 1))
cc <- ifelse(hdr %in% cols, NA, "NULL")   # "NULL" skips a column
base4 <- function() {
  d <- read.delim(f, quote = "", check.names = FALSE, colClasses = cc)
  d <- d[!(d$Decoy %in% "+") & d$PEP < 0.01, ]
  aggregate(Intensity ~ `Leading razor proteins`, data = d, na.action = na.pass,
            FUN = function(x) c(mean_int = mean(x, na.rm = TRUE), n = length(x)))
}
t_base4 <- replicate(5, system.time(base4())["elapsed"])
d <- fread(f, sep = "\t", quote = "", select = cols)
a <- d[!(Decoy %in% "+") & PEP < 0.01, .(mean_int = mean(Intensity, na.rm = TRUE), n = .N), by = `Leading razor proteins`]
cat(sprintf("rows %d | groups %d (base R, 4 columns: %d)\n", nrow(d), nrow(a), nrow(base4())))
cat(sprintf("R data.table fread (%d threads): median %.2f s\n", getDTthreads(), median(t_dt)))
cat(sprintf("R base read.delim + aggregate:   median %.2f s\n", median(t_base)))
cat(sprintf("R base read.delim, 4 columns:    median %.2f s\n", median(t_base4)))
