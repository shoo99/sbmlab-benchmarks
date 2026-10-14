## t-test vs limma vs DEqMS on the protein tables of three engines (PXD028735 HYE, 3 vs 3).
## Truth log2(B/A): human 0, yeast -1, E. coli +2. Proteins complete in both groups only.
## DEqMS count: min per-sample count + 1 (MaxQuant razor+unique peptides, FragPipe PSMs), DIA-NN N.Sequences + 1.
## Variant "median": per-sample median centring of log2 values before testing (no use of species labels).
suppressPackageStartupMessages({library(data.table); library(limma); library(DEqMS)})
## Contaminants: MaxQuant flags its own; FragPipe/Philosopher added none with a contam_ prefix and the DIA-NN FASTA had none,
## so protein groups containing any accession from MaxQuant 2.8.1 conf/contaminants.fasta are removed from FragPipe and DIA-NN.
con <- sub("^>", "", sub(" .*", "", grep("^>", readLines("contaminants_maxquant_2.8.1.fasta"), value = TRUE)))
A <- c("A_Alpha", "A_Beta", "A_Gamma"); B <- c("B_Alpha", "B_Beta", "B_Gamma")
sp1 <- function(x) { s <- unique(regmatches(x, gregexpr("_(HUMAN|YEAST|ECOLI)\\b", x))[[1]])
  s <- c(`_HUMAN` = "Human", `_YEAST` = "Yeast", `_ECOLI` = "Ecoli")[s]; if (length(s) == 1) unname(s) else NA_character_ }

## MaxQuant 2.8.1, LFQ, MBR on
pg <- fread("/home/sysoft/mq-hye/results_mbr_on/txt/proteinGroups.txt", sep = "\t", quote = "")
pg <- pg[!(Decoy %in% "+" | `Potential contaminant` %in% "+" | `Only identified by site` %in% "+")]
spmq <- function(h) { ox <- unique(regmatches(h, gregexpr("OX=[0-9]+", h))[[1]])
  s <- unique(c(`OX=9606` = "Human", `OX=559292` = "Yeast", `OX=83333` = "Ecoli")[ox]); s <- s[!is.na(s)]
  if (length(s) == 1) unname(s) else NA_character_ }
mq <- list(name = "MaxQuant 2.8.1 (LFQ)", species = vapply(pg$`Fasta headers`, spmq, ""),
           m = as.matrix(pg[, paste("LFQ intensity", c(A, B)), with = FALSE]),
           count = apply(as.matrix(pg[, paste("Razor + unique peptides", c(A, B)), with = FALSE]), 1, min))  # as in the DEqMS vignette
## FragPipe 24, MaxLFQ, MBR on
cp <- fread("/home/sysoft/fp-hye/work/combined_protein.tsv", sep = "\t", quote = "")
acc_fp <- function(p, ind) c(sub("^sp\\|([^|]+)\\|.*", "\\1", p), unlist(regmatches(ind, gregexpr("(?<=sp\\|)[^|]+", ind, perl = TRUE))))
fpc <- mapply(function(p, i) any(acc_fp(p, i) %in% con), cp$Protein, cp$`Indistinguishable Proteins`)
cat("FragPipe contaminant groups removed:", sum(fpc), "\n"); cp <- cp[!fpc]
fp <- list(name = "FragPipe 24 (MaxLFQ)", species = mapply(function(p, e) sp1(paste(p, e)), cp$Protein, cp$`Indistinguishable Proteins`),
           m = as.matrix(cp[, paste(c(A, B), "MaxLFQ Intensity"), with = FALSE]),
           count = apply(as.matrix(cp[, paste(c(A, B), "Spectral Count"), with = FALSE]), 1, min))  # PSMs: no per-sample peptide counts in combined_protein.tsv
## DIA-NN 2.3.2, pg_matrix (MaxLFQ), MBR on, SWATH
dn <- fread("/home/sysoft/blog-pipeline/eval/diann-hye/report.pg_matrix.tsv", sep = "\t")
dnc <- vapply(strsplit(dn$Protein.Group, ";"), function(a) any(a %in% con), TRUE)
cat("DIA-NN contaminant groups removed:", sum(dnc), "\n"); dn <- dn[!dnc]
runs <- grep("wiff$", names(dn), value = TRUE)
ord <- sapply(c("A_Sample_Alpha", "A_Sample_Beta", "A_Sample_Gamma", "B_Sample_Alpha", "B_Sample_Beta", "B_Sample_Gamma"), function(k) grep(k, runs, value = TRUE))
di <- list(name = "DIA-NN 2.3.2 (pg_matrix, SWATH)", species = vapply(dn$Protein.Names, sp1, ""),
           m = as.matrix(dn[, ..ord]), count = dn$N.Sequences)  # pg_matrix has no per-sample counts

run <- function(e, norm) {
  m <- e$m; m[m == 0] <- NA; m <- log2(m); colnames(m) <- c(A, B)
  keep <- !is.na(e$species) & rowSums(!is.na(m)) == 6
  m <- m[keep, ]; sp <- e$species[keep]; cnt <- e$count[keep]
  if (norm == "median") m <- sweep(m, 2, apply(m, 2, median))
  grp <- factor(rep(c("A", "B"), each = 3)); design <- model.matrix(~ grp)
  tt <- apply(m, 1, function(x) t.test(x[4:6], x[1:3])$p.value)
  fit <- eBayes(lmFit(m, design)); lim <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  fit$count <- cnt + 1; rownames(m) <- seq_len(nrow(m)); rownames(fit$coefficients) <- rownames(m)
  dq <- outputResult(spectraCounteBayes(fit), coef_col = 2); dq <- dq[order(as.integer(rownames(dq))), ]
  lfc <- rowMeans(m[, 4:6]) - rowMeans(m[, 1:3])
  adj <- list(t = p.adjust(tt, "BH"), limma = lim$adj.P.Val, DEqMS = dq$sca.adj.pva)
  rbindlist(lapply(names(adj), function(k) rbindlist(lapply(c(0, 0.5), function(fc) {
    call <- adj[[k]] < 0.05 & abs(lfc) >= fc
    data.table(engine = e$name, norm = norm, test = k, fc_min = fc, n_tested = nrow(m),
               n_human = sum(sp == "Human"), n_changed = sum(sp != "Human"),
               human_median_lfc = round(median(lfc[sp == "Human"]), 3),
               true_calls = sum(call & sp != "Human"),
               true_calls_right_sign = sum(call & ((sp == "Yeast" & lfc < 0) | (sp == "Ecoli" & lfc > 0))),
               false_calls = sum(call & sp == "Human"),
               false_up = sum(call & sp == "Human" & lfc > 0),
               FDP = round(100 * sum(call & sp == "Human") / max(sum(call), 1), 1))
  }))))
}
res <- rbindlist(lapply(list(mq, fp, di), function(e) rbind(run(e, "none"), run(e, "median"))))
fwrite(res, "stats_three_engines.csv")
print(res[, .(engine, norm, test, fc_min, n_tested, human_median_lfc, true_calls, n_changed, false_calls, FDP)], nrows = 100)
