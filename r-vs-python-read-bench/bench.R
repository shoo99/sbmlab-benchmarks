## Read MaxQuant evidence.txt (146 MB) + aggregate per protein; median of 5 runs.
suppressPackageStartupMessages(library(data.table))
f <- commandArgs(TRUE)[1]
t_dt <- replicate(5, system.time({
  d <- fread(f, sep = "\t", quote = "", select = c("Leading razor proteins", "PEP", "Intensity", "Decoy"))
  a <- d[!(Decoy %in% "+") & PEP < 0.01, .(mean_int = mean(Intensity, na.rm = TRUE), n = .N), by = `Leading razor proteins`]
})["elapsed"])
t_base <- replicate(5, system.time({
  d <- read.delim(f, quote = "", check.names = FALSE)
  d <- d[!(d$Decoy %in% "+") & d$PEP < 0.01, ]
  a <- aggregate(Intensity ~ `Leading razor proteins`, data = d, FUN = function(x) mean(x, na.rm = TRUE))
})["elapsed"])
d <- fread(f, sep = "\t", quote = "", select = c("Leading razor proteins", "PEP", "Intensity", "Decoy"))
a <- d[!(Decoy %in% "+") & PEP < 0.01, .(mean_int = mean(Intensity, na.rm = TRUE), n = .N), by = `Leading razor proteins`]
cat(sprintf("rows %d | groups %d\n", nrow(d), nrow(a)))
cat(sprintf("R data.table fread (%d threads): median %.2f s\n", getDTthreads(), median(t_dt)))
cat(sprintf("R base read.delim + aggregate:   median %.2f s\n", median(t_base)))
