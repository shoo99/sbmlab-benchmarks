## The real gene list of go-release-drift (issue #19): dexamethasone vs untreated airway smooth-muscle cells
## (Bioconductor 'airway', Himes et al. 2014, GSE52778), DESeq2 with the paired design ~ cell + condition,
## as in ../deseq2-edger-limma-airway (genes with >= 10 counts in >= 3 samples).
## Writes airway_de.tsv.gz: every tested gene with its log2 fold change and adjusted p-value. enrich.py takes
## padj < 0.05 and log2FC >= 1 (up) or <= -1 (down) as the two lists and all tested genes as their background.
## Usage: Rscript airway_de.R [data_dir]     (needs DESeq2; data/airway.RData from download.sh)
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2)})
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1] else "data"
load(file.path(data_dir, "airway.RData"))
counts <- assay(airway)
coldata <- as.data.frame(colData(airway))[, c("cell", "dex")]
coldata$condition <- factor(ifelse(coldata$dex == "trt", "treated", "control"), levels = c("control", "treated"))
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ cell + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "treated", "control"))
out <- data.frame(ensembl = rownames(res), log2fc = signif(res$log2FoldChange, 5), padj = signif(res$padj, 5))
write.table(out, gzfile("airway_de.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
up <- sum(out$padj < 0.05 & out$log2fc >= 1, na.rm = TRUE); down <- sum(out$padj < 0.05 & out$log2fc <= -1, na.rm = TRUE)
cat("genes tested", nrow(out), "| padj < 0.05 and |log2FC| >= 1:", up, "up,", down, "down\n")
cat("R", as.character(getRversion()), "| DESeq2", as.character(packageVersion("DESeq2")), "\n")
