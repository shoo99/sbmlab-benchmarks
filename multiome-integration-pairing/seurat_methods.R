## Unpaired RNA + ATAC integration scored with the true 10x Multiome pairing (issue #16): data preparation, the RNA
## cell types, and the two R methods (a naive baseline and Seurat's CCA label transfer / co-embedding).
##
## Data: 10x Genomics PBMC granulocyte-sorted 10k Multiome (Cell Ranger ARC 2.0.0): RNA and ATAC from the same nuclei,
## so the true match of every ATAC profile is the RNA profile with the same barcode. The methods never see this pairing;
## it is used only to score them.
## Cells: RNA counts 1,000-25,000, at least 500 genes, mitochondrial share < 20 %, ATAC peak counts 1,000-100,000.
## RNA cell types (the label truth): log-normalise, 2,000 variable genes, 30 PCs, Louvain clusters (resolution 0.8);
##   each cluster gets the cell type whose marker set has the highest mean scaled expression (markers below); clusters
##   whose best score is below 0.15 are left out as unassigned.
## Gene activity: ATAC peak counts in each gene's body plus 2 kb upstream (strand-aware; gene spans from the iGenomes
##   UCSC hg38 genes.gtf, matched to the 10x gene symbols).
## Methods
##   naive baseline       RNA and gene activity on the shared variable genes, each z-scored within its modality,
##                        one PCA (30 PCs) of both; labels by majority of the 15 nearest RNA cells
##   Seurat CCA           FindTransferAnchors(reduction = "cca") from RNA to gene activity, TransferData for labels
##                        (weights from the ATAC LSI, TF-IDF + SVD of peaks, dims 2-30); co-embedding as in Seurat's
##                        ATAC vignette: RNA imputed for ATAC cells from the anchors, merged with RNA, one PCA (30 PCs)
## Writes work/ (inputs for run_glue.py) and results/<method>.tsv: one row per cell with its true type, the transferred
## type, and FOSCTTM in both directions (share of the other modality's cells closer than the true match).
## Usage: Rscript seurat_methods.R [data_dir] > seurat_methods.out     (seed 1 throughout)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(hdf5r); library(data.table) })
options(width = 140)
set.seed(1)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
args <- commandArgs(TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path(script_dir, "data")
work <- file.path(script_dir, "work"); dir.create(work, showWarnings = FALSE)
res_dir <- file.path(script_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }

## ---- read the h5 with feature coordinates
h5 <- H5File$new(file.path(data_dir, "pbmc10k_multiome.h5"), "r")
mx <- sparseMatrix(i = h5[["matrix/indices"]][] + 1, p = h5[["matrix/indptr"]][], x = as.numeric(h5[["matrix/data"]][]),
                   dims = h5[["matrix/shape"]][], repr = "C")
ftype <- h5[["matrix/features/feature_type"]][]; fname <- h5[["matrix/features/name"]][]
fint <- h5[["matrix/features/interval"]][]; barcodes <- h5[["matrix/barcodes"]][]
h5$close_all()
colnames(mx) <- barcodes
rna <- mx[ftype == "Gene Expression", ]; rownames(rna) <- make.unique(fname[ftype == "Gene Expression"])
tss <- data.frame(gene = rownames(rna), interval = fint[ftype == "Gene Expression"])   # 10x TSS region of each gene
peaks <- mx[ftype == "Peaks", ]; rownames(peaks) <- fint[ftype == "Peaks"]
cat(sprintf("== read: %d barcodes, %d genes, %d peaks\n", ncol(mx), nrow(rna), nrow(peaks)))
rm(mx); invisible(gc())

## ---- cell filter
n_rna <- Matrix::colSums(rna); g_rna <- Matrix::colSums(rna > 0); n_atac <- Matrix::colSums(peaks)
mt <- Matrix::colSums(rna[grepl("^MT-", rownames(rna)), ]) / pmax(n_rna, 1)
keep <- n_rna >= 1000 & n_rna <= 25000 & g_rna >= 500 & mt < 0.2 & n_atac >= 1000 & n_atac <= 1e5
rna <- rna[, keep]; peaks <- peaks[, keep]
cat(sprintf("   kept %d cells (RNA 1,000-25,000 counts, >= 500 genes, mito < 20 %%, ATAC 1,000-100,000)\n", ncol(rna)))

## ---- RNA cell types
ro <- quiet(CreateSeuratObject(rna, min.cells = 3))
ro <- quiet(NormalizeData(ro)); ro <- quiet(FindVariableFeatures(ro, nfeatures = 2000))
MARKER_GENES <- c("CD14", "LYZ", "S100A8", "VCAN", "FCGR3A", "MS4A7", "CDKN1C", "FCER1A", "CD1C", "CLEC10A", "LILRA4",
                  "IL3RA", "CLEC4C", "MS4A1", "CD79A", "BANK1", "PAX5", "CD3E", "IL7R", "CD4", "CCR7", "CD8A", "CD8B",
                  "NKG7", "GNLY", "KLRD1", "PRF1", "PPBP", "PF4")
ro <- quiet(ScaleData(ro, features = union(VariableFeatures(ro), intersect(MARKER_GENES, rownames(ro)))))
ro <- quiet(RunPCA(ro, npcs = 30, seed.use = 1))
ro <- quiet(FindNeighbors(ro, dims = 1:30)); ro <- quiet(FindClusters(ro, resolution = 0.8, random.seed = 1))
MARKERS <- list("CD14 Mono" = c("CD14", "LYZ", "S100A8", "VCAN"), "CD16 Mono" = c("FCGR3A", "MS4A7", "CDKN1C"),
                "cDC" = c("FCER1A", "CD1C", "CLEC10A"), "pDC" = c("LILRA4", "IL3RA", "CLEC4C"),
                "B" = c("MS4A1", "CD79A", "BANK1", "PAX5"), "CD4 T" = c("CD3E", "IL7R", "CD4", "CCR7"),
                "CD8 T" = c("CD3E", "CD8A", "CD8B"), "NK" = c("NKG7", "GNLY", "KLRD1", "PRF1"),
                "Platelet" = c("PPBP", "PF4"))
sc <- GetAssayData(ro, layer = "scale.data")
cl <- Idents(ro)
score <- sapply(MARKERS, function(g) tapply(colMeans(sc[intersect(g, rownames(sc)), , drop = FALSE]), cl, mean))
## T cells: CD4 vs CD8 by the CD8 genes, NK vs T by CD3E (both sets share CD3E or NKG7-like genes)
cell_type_of_cluster <- ifelse(apply(score, 1, max) >= 0.15, colnames(score)[apply(score, 1, which.max)], "unassigned")
names(cell_type_of_cluster) <- rownames(score)
cat("\n== RNA clusters: mean scaled marker expression and the assigned type (unassigned: best score < 0.15)\n")
print(data.frame(round(score, 2), cells = as.numeric(table(cl)), type = cell_type_of_cluster, check.names = FALSE))
ct_all <- unname(cell_type_of_cluster[as.character(cl)])
assigned <- ct_all != "unassigned"
rna <- rna[, assigned]; peaks <- peaks[, assigned]; ro <- ro[, assigned]
ct <- factor(ct_all[assigned])
cat(sprintf("\n   %d cells in unassigned clusters dropped; %d cells left\n", sum(!assigned), ncol(rna)))
cat("\n   cells per type:\n"); print(table(ct))

## ---- gene activity: peaks in gene body + 2 kb upstream
## gene bodies: the transcripts (iGenomes UCSC hg38 genes.gtf) on the chromosome of the gene's 10x TSS region that
## overlap it; the GTF repeats some gene names at distant loci, so a plain min-max per name would span megabases
gtf <- fread(file.path(data_dir, "hg38_genes.gtf"), sep = "\t", header = FALSE, quote = "", select = c(1, 4, 5, 7, 9),
             col.names = c("chr", "start", "end", "strand", "attr"))
gtf[, gene := sub('.*gene_name "([^"]+)".*', "\\1", attr)][, tx := sub('.*transcript_id "([^"]+)".*', "\\1", attr)][, attr := NULL]
tx <- gtf[, .(start = min(start), end = max(end), strand = strand[1]), by = .(gene, tx, chr)]
rm(gtf); invisible(gc())
tp <- tss[grepl("^chr", tss$interval), ]
tp <- data.table(gene = tp$gene, tchr = sub(":.*", "", tp$interval), ts = as.integer(sub(".*:(\\d+)-.*", "\\1", tp$interval)),
                 te = as.integer(sub(".*-", "", tp$interval)))
tx <- merge(tx, tp, by = "gene")[chr == tchr & start <= te & end >= ts]
gs <- as.data.frame(tx[, .(chr = chr[1], start = min(start), end = max(end), strand = strand[1]), by = gene])
gs <- gs[gs$gene %in% rownames(rna) & grepl("^chr([0-9]+|X|Y)$", gs$chr), ]
gs$ws <- ifelse(gs$strand == "+", gs$start - 2000, gs$start); gs$we <- ifelse(gs$strand == "+", gs$end, gs$end + 2000)
pk <- do.call(rbind, strsplit(rownames(peaks), "[:-]"))
pk <- data.frame(chr = pk[, 1], start = as.integer(pk[, 2]), end = as.integer(pk[, 3]), i = seq_len(nrow(pk)))
ii <- integer(0); jj <- integer(0)
for (chr in unique(gs$chr)) {
  p <- pk[pk$chr == chr, ]; p <- p[order(p$start), ]; g <- which(gs$chr == chr)
  if (!nrow(p)) next
  lo <- findInterval(gs$ws[g] - 1, p$end) + 1            # first peak ending at or after the window start
  hi <- findInterval(gs$we[g], p$start)                  # last peak starting at or before the window end
  for (k in seq_along(g)) if (hi[k] >= lo[k]) { ii <- c(ii, rep(g[k], hi[k] - lo[k] + 1)); jj <- c(jj, p$i[lo[k]:hi[k]]) }
}
G2P <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(nrow(gs), nrow(peaks)), dimnames = list(gs$gene, NULL))
ga <- G2P %*% peaks
cat(sprintf("\n== gene activity: %d genes with coordinates, %d gene-peak links, median %d peak counts per cell in gene windows\n",
            nrow(gs), length(ii), as.integer(median(Matrix::colSums(ga)))))

## ATAC LSI (TF-IDF of peak counts, then 30 SVD components), computed in place to keep memory down; the peak matrix
## is not needed after this
cs <- Matrix::colSums(peaks); df_ <- tabulate(peaks@i + 1, nrow(peaks)); idf <- log1p(ncol(peaks) / pmax(df_, 1))
peaks@x <- log1p(peaks@x / cs[rep.int(seq_len(ncol(peaks)), diff(peaks@p))] * idf[peaks@i + 1] * 1e4)
lsi <- irlba::irlba(peaks, nv = 30)
lsi_emb <- lsi$v %*% diag(lsi$d); rownames(lsi_emb) <- colnames(peaks); colnames(lsi_emb) <- paste0("LSI_", 1:30)
rm(peaks, G2P, lsi); invisible(gc())

## ---- inputs for run_glue.py: the cells and their types, and gene coordinates (it reads the counts from the h5)
utils::write.table(gs[, c("gene", "chr", "start", "end", "strand")], file.path(work, "gene_coords.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(barcode = colnames(rna), cell_type = ct), file.path(work, "cells.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

## ---- scoring helpers
foscttm <- function(A, R) {   # rows: cells in the same order (true pairs); returns per-cell FOSCTTM both ways
  n <- nrow(A); fa <- fr <- numeric(n); na <- rowSums(A^2); nr <- rowSums(R^2)
  for (s in split(seq_len(n), ceiling(seq_len(n) / 1000))) {
    self <- cbind(seq_along(s), s)       # the true pair's distance from the same matrix, so rounding cannot count it as closer
    d <- outer(na[s], nr, "+") - 2 * A[s, , drop = FALSE] %*% t(R)                   # ATAC cells s vs all RNA
    fa[s] <- rowSums(d < d[self]) / (n - 1)
    d2 <- outer(nr[s], na, "+") - 2 * R[s, , drop = FALSE] %*% t(A)                  # RNA cells s vs all ATAC
    fr[s] <- rowSums(d2 < d2[self]) / (n - 1)
  }
  list(atac = fa, rna = fr)
}
knn_labels <- function(A, R, labels, k = 15) {
  out <- character(nrow(A)); nr <- rowSums(R^2)
  for (s in split(seq_len(nrow(A)), ceiling(seq_len(nrow(A)) / 1000))) {
    d <- outer(rowSums(A[s, , drop = FALSE]^2), nr, "+") - 2 * A[s, , drop = FALSE] %*% t(R)
    out[s] <- apply(d, 1, function(x) { t_ <- table(labels[order(x)[1:k]]); names(t_)[which.max(t_)] })
  }
  out
}
save_result <- function(method, pred, f) {
  d <- data.frame(barcode = colnames(rna), true_type = ct, predicted_type = pred,
                  foscttm_atac = round(f$atac, 6), foscttm_rna = round(f$rna, 6))
  utils::write.table(d, file.path(res_dir, paste0(method, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("   %s: FOSCTTM %.4f (ATAC->RNA %.4f, RNA->ATAC %.4f), label accuracy %.3f\n", method,
              mean(c(f$atac, f$rna)), mean(f$atac), mean(f$rna), mean(pred == as.character(ct))))
}

## ---- naive baseline
cat("\n== methods\n")
genes <- intersect(VariableFeatures(ro), rownames(ga))
ln <- function(m) log1p(t(t(m) / pmax(Matrix::colSums(m), 1)) * 1e4)
zr <- t(scale(t(as.matrix(ln(rna[genes, ]))))); za <- t(scale(t(as.matrix(ln(ga[genes, ])))))
zr[is.na(zr)] <- 0; za[is.na(za)] <- 0
pc <- irlba::prcomp_irlba(t(cbind(zr, za)), n = 30, center = TRUE, scale. = FALSE)$x
n <- ncol(rna); R <- pc[1:n, ]; A <- pc[n + seq_len(n), ]
save_result("naive_baseline", knn_labels(A, R, as.character(ct)), foscttm(A, R))

## ---- Seurat CCA
go <- quiet(CreateSeuratObject(ga[, colnames(rna)], assay = "RNA"))
go <- quiet(NormalizeData(go)); go <- quiet(ScaleData(go, features = genes))
go[["lsi"]] <- CreateDimReducObject(lsi_emb, key = "LSI_", assay = "RNA")
anchors <- quiet(FindTransferAnchors(reference = ro, query = go, features = genes, reference.assay = "RNA",
                                     query.assay = "RNA", reduction = "cca"))
pred <- quiet(TransferData(anchors, refdata = as.character(ct), weight.reduction = go[["lsi"]], dims = 2:30))
imp <- quiet(TransferData(anchors, refdata = GetAssayData(ro, layer = "data")[genes, ], weight.reduction = go[["lsi"]],
                          dims = 2:30))
imputed <- as.matrix(GetAssayData(imp, layer = "data"))
comb <- cbind(as.matrix(GetAssayData(ro, layer = "data")[genes, ]), imputed)
colnames(comb) <- c(paste0("rna_", colnames(rna)), paste0("atac_", colnames(rna)))
co <- t(scale(t(comb))); co[is.na(co)] <- 0
pc2 <- irlba::prcomp_irlba(t(co), n = 30, center = TRUE, scale. = FALSE)$x
save_result("seurat_cca", as.character(pred$predicted.id), foscttm(pc2[n + seq_len(n), ], pc2[1:n, ]))

cat("\n", R.version.string, " | Seurat ", format(packageVersion("Seurat")), " | Matrix ", format(packageVersion("Matrix")),
    " | irlba ", format(packageVersion("irlba")), "\n", sep = "")
