## The sbmlab BioMart tutorial's R code, run for Ensembl 111 (the post's pin) and the current release.
suppressPackageStartupMessages(library(biomaRt))
cat("biomaRt", as.character(packageVersion("biomaRt")), "\n")
run <- function(v) {
  ## useEnsembl(version = 111) failed in biomaRt 2.58.2 ("invalid 'file' argument" while listing archives),
  ## so point at the archive host directly. Release 111 = January 2024 archive.
  host <- if (is.null(v)) "https://www.ensembl.org" else "https://jan2024.archive.ensembl.org"
  pig <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "sscrofa_gene_ensembl", host = host)
  orthos <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene",
                                 "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type",
                                 "mmusculus_homolog_perc_id"), mart = pig)
  orthos <- orthos[orthos$mmusculus_homolog_ensembl_gene != "", ]
  cat(sprintf("\nEnsembl %s: rows with a mouse homolog %d\n", ifelse(is.null(v), "current", v), nrow(orthos)))
  print(table(orthos$mmusculus_homolog_orthology_type))
  one <- subset(orthos, mmusculus_homolog_orthology_type == "ortholog_one2one")
  cat("one2one pairs:", nrow(one), " unique pig genes:", length(unique(one$ensembl_gene_id)),
      " with a pig gene name:", sum(one$external_gene_name != ""), "\n")
  write.table(one, sprintf("one2one_%s.tsv", ifelse(is.null(v), "current", v)), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(one)
}
a <- run(111); b <- run(NULL)
k <- function(x) paste(x$ensembl_gene_id, x$mmusculus_homolog_ensembl_gene)
cat("\nshared pairs 111 vs current:", length(intersect(k(a), k(b))), "\n")
