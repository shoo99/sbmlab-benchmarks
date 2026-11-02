## Query BioMart's martservice directly from R (no biomaRt), release 111 archive.
attrs <- c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene",
           "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type", "mmusculus_homolog_perc_id")
xml <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
  '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" completionStamp="1">',
  '<Dataset name="sscrofa_gene_ensembl" interface="default">',
  paste0('<Attribute name="', attrs, '"/>', collapse = ""), '</Dataset></Query>')
url <- paste0("https://jan2024.archive.ensembl.org/biomart/martservice?query=", URLencode(xml, reserved = TRUE))
## added 2026-11-02 (not in the post): BioMart can cut a result short or answer "Query ERROR ..." with HTTP 200;
## completionStamp="1" ends a complete result with "[success]", so check for it before parsing
txt <- readLines(url, warn = FALSE)
if (!length(txt) || tail(txt, 1) != "[success]") stop("BioMart result incomplete or an error: ", head(txt, 1))
orthos <- read.delim(text = head(txt, -1), header = FALSE, col.names = attrs, colClasses = "character")
orthos <- orthos[orthos$mmusculus_homolog_ensembl_gene != "", ]
one2one <- subset(orthos, mmusculus_homolog_orthology_type == "ortholog_one2one")
cat("rows with mouse homolog", nrow(orthos), " one2one", nrow(one2one), "\n")
