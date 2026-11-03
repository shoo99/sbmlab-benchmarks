## Gene-level count matrices (genes x samples) from each tool's output in <work dir> (written by run.sh).
## Salmon NumReads and kallisto est_counts are summed over the transcripts of each gene (transcript -> gene from the
## GTF; ERCC sequences are their own gene); RSEM gives expected_count per gene; featureCounts counts per gene.
## Usage: Rscript gene_matrix.R <work dir> [out dir = .]
## Writes genes_<tool>.tsv.gz (gene_id + one column per sample), and copies versions.tsv and timings.tsv.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
W <- a[1]; OUT <- if (length(a) > 1) a[2] else "."
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
samples <- readLines(file.path(W, "samples.txt"))

gtf <- fread(file.path(W, "annotation_ercc.gtf"), sep = "\t", header = FALSE, quote = "")[V3 == "exon"]
attr_of <- function(x, key) sub(sprintf('.*%s "([^"]+)".*', key), "\\1", x)
tx2gene <- unique(data.table(tx = attr_of(gtf$V9, "transcript_id"), gene = attr_of(gtf$V9, "gene_id")))
gene_of <- function(tx) {
  g <- tx2gene$gene[match(tx, tx2gene$tx)]
  g[grepl("^ERCC-", tx)] <- tx[grepl("^ERCC-", tx)]            # ERCC FASTA names are the gene IDs
  if (anyNA(g)) stop(sum(is.na(g)), " transcripts without a gene, e.g. ", tx[is.na(g)][1])
  g
}
per_sample <- function(read_one) {
  m <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), lapply(samples, function(s) {
    d <- read_one(s)
    setnames(d[, .(count = sum(count)), by = gene_id], "count", s)
  }))
  m[is.na(m)] <- 0
  m
}

mats <- list(
  salmon = per_sample(function(s) {
    q <- fread(file.path(W, "salmon", s, "quant.sf"))
    data.table(gene_id = gene_of(q$Name), count = q$NumReads)
  }),
  kallisto = per_sample(function(s) {
    q <- fread(file.path(W, "kallisto", s, "abundance.tsv"))
    data.table(gene_id = gene_of(q$target_id), count = q$est_counts)
  }),
  rsem = per_sample(function(s) {
    q <- fread(file.path(W, "rsem", paste0(s, ".genes.results")))
    data.table(gene_id = q$gene_id, count = q$expected_count)
  }),
  featurecounts = {
    fc <- fread(file.path(W, "featurecounts", "counts.txt"), skip = 1)
    m <- fc[, c(1, 7:ncol(fc)), with = FALSE]
    setnames(m, c("gene_id", samples))                          # run.sh passes the BAMs in samples.txt order
    m
  })

for (tool in names(mats)) {
  m <- mats[[tool]]
  setcolorder(m, c("gene_id", samples))
  setorder(m, gene_id)
  num <- samples
  m[, (num) := lapply(.SD, function(v) round(v, 2)), .SDcols = num]
  fwrite(m, file.path(OUT, sprintf("genes_%s.tsv.gz", tool)), sep = "\t")
  cat(sprintf("%-14s %6d genes (%d ERCC) x %d samples\n", tool, nrow(m), sum(grepl("^ERCC-", m$gene_id)), length(samples)))
}
for (f in c("versions.tsv", "timings.tsv")) file.copy(file.path(W, f), file.path(OUT, f), overwrite = TRUE)
