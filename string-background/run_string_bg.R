## STRING enrichment background test (2026-10-13).
## Null: random 300 proteins drawn from the proteins actually detected in a DIA-NN run (HYE human part).
## Planted: 15 detected cholesterol-biosynthesis proteins + 285 random detected proteins.
## Each list is enriched twice: STRING's default background (whole genome) vs the detected proteins.
suppressPackageStartupMessages({library(data.table); library(STRINGdb); library(msigdbr)})
CACHE <- "string_cache"

## 1. detected human proteins -> STRING ids
pg <- fread("../diann-hye/report.pg_matrix.tsv")
det <- pg[grepl("_HUMAN", Protein.Names) & !grepl("_YEAST|_ECOLI", Protein.Names),
          .(acc = sub(";.*", "", Protein.Group), gene = sub(";.*", "", Genes))]
sdb <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = CACHE)
det <- as.data.table(sdb$map(as.data.frame(det), "acc", removeUnmappedRows = TRUE))
det <- det[!duplicated(STRING_id)]
cat("detected human proteins mapped to STRING:", nrow(det), "\n")

## 2. enrichment helper: count terms with FDR < 0.05 per category
enrich <- function(ids, background = NULL) {
  s <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = CACHE)
  if (!is.null(background)) s$set_background(background)
  e <- as.data.table(s$get_enrichment(ids))
  Sys.sleep(1)
  if (!nrow(e)) return(e)
  e[fdr < 0.05]
}
count <- function(e) {
  cats <- c("Process", "Component", "Function", "KEGG")
  setNames(vapply(cats, function(k) if (nrow(e)) sum(e$category == k) else 0L, 1L), cats)
}

## 3. null lists
res <- list(); terms <- list()
for (seed in 1:20) {
  set.seed(seed); ids <- sample(det$STRING_id, 300)
  g <- enrich(ids); d <- enrich(ids, det$STRING_id)
  res[[length(res) + 1]] <- data.table(list = "random", seed, background = c("genome", "detected"),
                                       rbind(count(g), count(d)))
  terms[[seed]] <- if (nrow(g)) g[, .(seed, category, term, description, fdr)] else NULL
  cat("seed", seed, "genome:", nrow(g), " detected:", nrow(d), "\n")
}

## 4. planted pathway
ms <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
chol <- unique(ms$gene_symbol[ms$gs_name == "REACTOME_CHOLESTEROL_BIOSYNTHESIS"])
chol_ids <- det[gene %in% chol, STRING_id]
cat("cholesterol-biosynthesis proteins detected:", length(chol_ids), "\n")
planted <- list()
for (seed in 1:5) {
  set.seed(100 + seed)
  ids <- c(sample(chol_ids, 15), sample(setdiff(det$STRING_id, chol_ids), 285))
  g <- enrich(ids); d <- enrich(ids, det$STRING_id)
  hit <- function(e) if (nrow(e)) any(grepl("holesterol|terol biosynth", e$description)) else FALSE
  res[[length(res) + 1]] <- data.table(list = "planted", seed, background = c("genome", "detected"),
                                       rbind(count(g), count(d)))
  planted[[seed]] <- data.table(seed, found_genome = hit(g), found_detected = hit(d),
                                n_genome = nrow(g), n_detected = nrow(d))
}

res <- rbindlist(res); fwrite(res, "counts.csv"); print(res)
fwrite(rbindlist(terms), "null_terms_genome_background.csv")
print(rbindlist(planted)); fwrite(rbindlist(planted), "planted.csv")
