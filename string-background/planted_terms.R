## Which terms come back for the planted list (seed 101), per background
suppressPackageStartupMessages({library(data.table); library(STRINGdb); library(msigdbr)})
CACHE <- "string_cache"
pg <- fread("../diann-hye/report.pg_matrix.tsv")
det <- pg[grepl("_HUMAN", Protein.Names) & !grepl("_YEAST|_ECOLI", Protein.Names), .(acc = sub(";.*", "", Protein.Group), gene = sub(";.*", "", Genes))]
s <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = CACHE)
det <- as.data.table(s$map(as.data.frame(det), "acc", removeUnmappedRows = TRUE)); det <- det[!duplicated(STRING_id)]
ms <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
chol_ids <- det[gene %in% unique(ms$gene_symbol[ms$gs_name == "REACTOME_CHOLESTEROL_BIOSYNTHESIS"]), STRING_id]
set.seed(101); ids <- c(sample(chol_ids, 15), sample(setdiff(det$STRING_id, chol_ids), 285))
g <- as.data.table(s$get_enrichment(ids)); invisible(s$set_background(det$STRING_id)); d <- as.data.table(s$get_enrichment(ids))
k <- c("Process", "Component", "Function", "KEGG")
cat("detected background, GO/KEGG FDR<0.05:\n"); print(d[fdr < 0.05 & category %in% k, .(category, description, number_of_genes, fdr)])
cat("genome background, GO/KEGG FDR<0.05:", nrow(g[fdr < 0.05 & category %in% k]), "; of which mention sterol/cholesterol/lipid:",
    nrow(g[fdr < 0.05 & category %in% k & grepl("terol|lipid|isoprenoid|steroid", description, ignore.case = TRUE)]), "\n")
gg <- g[fdr < 0.05 & category %in% k][order(fdr)]; gg[, rank := .I]
cat("genome background: rank of first sterol/cholesterol term:", gg[grepl("terol", description), min(rank)], "of", nrow(gg), "\n")
print(gg[1:10, .(rank, category, description, fdr)])
