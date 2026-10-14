library(limma)
library(DEqMS)
library(readr)

# Load MaxQuant proteinGroups.txt (or DIA-NN pg_matrix.tsv with peptide count joined)
pg <- read_tsv("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", show_col_types = FALSE)

# Standard filters. MaxQuant 2.8 renamed "Reverse" to "Decoy"; on a 2.8 file the old
# pg$Reverse != "+" filter silently returns zero rows. This version works on both:
flag <- function(col) if (col %in% names(pg)) pg[[col]] %in% "+" else FALSE
pg <- pg[!(flag("Reverse") | flag("Decoy") | flag("Potential contaminant") | flag("Only identified by site")), ]

# Identify intensity columns (LFQ or iBAQ)
lfq_cols <- grep("^LFQ intensity ", colnames(pg), value = TRUE)
mat <- as.matrix(pg[, lfq_cols])
rownames(mat) <- pg$`Majority protein IDs`

# Replace zeros with NA, log2 transform
mat[mat == 0] <- NA
mat <- log2(mat)

# Group assignments
group <- factor(c("Ctrl","Ctrl","Ctrl","Treat","Treat","Treat"))  # HYE: A = Ctrl, B = Treat
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(Treat - Ctrl, levels = design)

fit <- lmFit(mat, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)
limma_results <- topTable(fit, number = Inf, adjust.method = "BH")

fit <- lmFit(mat, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# DEqMS extension: peptide count per protein, aligned to the rows of the fit.
# As in the DEqMS vignette: the minimum "Razor + unique peptides" across samples, plus 1.
pep <- as.matrix(pg[, grep("^Razor \\+ unique peptides ", colnames(pg))])
rownames(pep) <- pg$`Majority protein IDs`
fit$count <- apply(pep, 1, min)[rownames(fit$coefficients)] + 1
deq <- spectraCounteBayes(fit)
deqms_results <- outputResult(deq, coef_col = 1)

cat("rows after filters", nrow(pg), "\nlimma rows", nrow(limma_results), " sig", sum(limma_results$adj.P.Val<0.05, na.rm=TRUE), "\nDEqMS rows", nrow(deqms_results), " sig", sum(deqms_results$sca.adj.pva<0.05, na.rm=TRUE), "\n")
