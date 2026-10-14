library(data.table); library(STRINGdb)

# 1. Proteins actually detected in the experiment (here: human protein groups in a DIA-NN pg_matrix)
pg  <- fread("report.pg_matrix.tsv")
det <- pg[grepl("_HUMAN", Protein.Names) & !grepl("_YEAST|_ECOLI", Protein.Names),
          .(acc = sub(";.*", "", Protein.Group))]

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400,
                          input_directory = "string_cache")
det <- as.data.table(string_db$map(as.data.frame(det), "acc", removeUnmappedRows = TRUE))
background <- unique(det$STRING_id)

# 2. A list with no biology in it: 300 detected proteins picked at random
set.seed(1)
hits <- sample(background, 300)

# 3. Same list, two backgrounds
genome_bg <- string_db$get_enrichment(hits)                 # STRING default: whole genome

invisible(string_db$set_background(background))             # detected proteins only
detected_bg <- string_db$get_enrichment(hits)

go_kegg <- c("Process", "Component", "Function", "KEGG")
sig <- function(e) subset(e, fdr < 0.05 & category %in% go_kegg)
c(genome = nrow(sig(genome_bg)), detected = nrow(sig(detected_bg)))
head(sig(genome_bg)[order(sig(genome_bg)$fdr), c("category", "description", "fdr")], 8)
