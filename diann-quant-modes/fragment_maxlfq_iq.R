## Fragment-level MaxLFQ with iq (Pham et al. 2020) on DIA-NN --export-quant fragment apex quantities (normalised).
suppressPackageStartupMessages({library(data.table); library(iq)})
a <- commandArgs(TRUE); if (length(a)) setwd(a[1])   # e.g. clean/nombr (default: this folder = first reruns)
d <- fread("fragments_long.tsv")
setnames(d, c("sample", "Protein.Group", "Protein.Names", "ion", "q"), c("sample", "protein", "names", "ion", "q"))
nm <- unique(d[, .(protein, names)])
pre <- iq::preprocess(as.data.frame(d[, .(protein, sample, ion, q)]), primary_id = "protein", sample_id = "sample",
                      secondary_id = "ion", intensity_col = "q", median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL)
pl <- iq::create_protein_list(pre)
res <- iq::create_protein_table(pl)
est <- as.data.table(2^res$estimate, keep.rownames = "Protein.Group")
est <- merge(nm[, .(Protein.Group = protein, Protein.Names = names)], est, by = "Protein.Group")
fwrite(est, "pg_fragment_maxlfq_iq.tsv", sep = "\t")
cat("proteins:", nrow(est), "\n")
