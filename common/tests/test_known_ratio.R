## Tests for common/known_ratio.R (base R only): a hand-worked toy matrix, then re-creating the bundled
## diann-hye/diann_scores.csv and diann-2.7-hye/diann_scores_270.csv with the shared functions.
## Run from anywhere: Rscript common/tests/test_known_ratio.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(root, "common", "known_ratio.R"))
fails <- 0
check <- function(desc, ok) {
  if (isTRUE(ok)) cat("ok  ", desc, "\n") else { cat("FAIL", desc, "\n"); fails <<- fails + 1 }
}

## ---- toy matrix: columns A1..A3, B1..B3 (linear); log2 values chosen so the expected results can be worked out
L <- rbind(c(10, 10, 10, 10, 10, 10),       # H1  Human  fc 0
           c(10, 10, 10, 10.4, 10.4, 10.4), # H2  Human  fc 0.4
           c(12, 12, 12, 11, 11, 11),       # Y1  Yeast  fc -1
           c(8, 8, NA, 10.6, 10.6, 10.6),   # E1  E. coli, 2 of 3 in A -> scored, fc 2.6
           c(8, NA, NA, 10, 10, 10),        # E2  E. coli, 1 of 3 in A -> not scored
           c(5, 5, 5, 5, 5, 5),             # two species in one group -> no species
           c(5, 5, 5, 5, 5, 5),             # "_HUMANX" is not "_HUMAN"
           c(10, 10, 10, 11, 11, 11))       # Y2  Yeast  fc +1 (wrong direction if called)
m <- 2^L
m[is.na(m)] <- 0                            # zeros must become missing
colnames(m) <- c("A1", "A2", "A3", "B1", "B2", "B3")
nm <- c("H1_HUMAN", "H2_HUMAN", "Y1_YEAST", "E1_ECOLI", "E2_ECOLI", "X_HUMAN;Y_YEAST", "Z_HUMANX", "Y2_YEAST")

sp <- species_from_names(nm)
check("species_from_names", identical(sp, c("Human", "Human", "Yeast", "Ecoli", "Ecoli", NA, NA, "Yeast")))
check("species_from_taxa", identical(species_from_taxa(c("OX=9606", "OX=9606 OX=559292", "OX=10090", "OX=10090;OX=83333", "")),
                                     c("Human", NA, NA, "Ecoli", NA)))

r <- log2_ratio(m, a = 1:3, b = 4:6)
check("log2_ratio ok (2 of 3 per group)", identical(r$ok, c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)))
check("log2_ratio fc", isTRUE(all.equal(r$fc, c(0, 0.4, -1, 2.6, 2, 0, 0, 1))))
check("log2_ratio by column name", identical(log2_ratio(m, c("A1", "A2", "A3"), c("B1", "B2", "B3")), r))

g <- group_log2_means(m, list(A = 1:3, B = 4:6))
check("group_log2_means", isTRUE(all.equal(unname(g[4:5, "A"]), c(8, NA))) && isTRUE(all.equal(unname(g[2, "B"]), 10.4)))

s <- score_groups(r$fc, r$ok, sp, HYE_TRUTH, centre = "Human")
## Human: fc 0, 0.4 -> median 0.2, IQR 0.2 (type 7), both within 0.5; centre shift = 0.2
## Yeast: Y1 -1 and Y2 +1 -> median 0, within 1 of 2 (Y2 is 2 away); centred -1.2 and 0.8 -> 1 of 2
## E. coli: only E1 (2.6): uncentred 0.6 away -> 0; centred 2.4 -> 0.4 away -> 1
check("score_groups n", identical(s$n, c(2L, 2L, 1L)))
check("score_groups median", isTRUE(all.equal(s$median, c(0.2, 0, 2.6))))
check("score_groups IQR", isTRUE(all.equal(s$IQR, c(0.2, 1, 0))))
check("score_groups within", isTRUE(all.equal(s$within, c(1, 0.5, 0))))
check("score_groups within_centred", isTRUE(all.equal(s$within_centred, c(1, 0.5, 1))))
check("score_groups centre median", isTRUE(all.equal(attr(s, "centre_median"), 0.2)))
check("score_groups without centre has no within_centred", !("within_centred" %in% names(score_groups(r$fc, r$ok, sp, HYE_TRUTH))))

called <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
o <- call_outcome(called, r$fc, sp, HYE_TRUTH)
check("call_outcome", identical(o, c("false", "not_called", "right", "right", "right", "no_truth", "not_called", "wrong_direction")))
cs <- call_summary(o)
check("call_summary", cs$calls == 6 && cs$right == 3 && cs$wrong_direction == 1 && cs$false == 1 && isTRUE(all.equal(cs$fdp, 0.4)))

tmp <- tempfile(fileext = ".csv")
write_scores(scores_rows("f", "d", "t", "1.0", "default", c("Human", "Yeast"), "within_0.5", c(1, 0.5), c(2, 2)), tmp)
check("write_scores header", identical(names(read.csv(tmp)), SCORES_COLUMNS))

## ---- re-create the bundled DIA-NN score tables (same rules as diann-hye/score_diann.R)
recreate <- function(dir, files, labels, ref) {
  res <- do.call(rbind, Map(function(f, label) {
    x <- utils::read.delim(file.path(root, dir, f), check.names = FALSE, quote = "", comment.char = "")
    rc <- grep("Condition_", names(x), value = TRUE)
    A <- grep("Condition_A", rc, value = TRUE); B <- grep("Condition_B", rc, value = TRUE)
    r <- log2_ratio(x[, c(A, B)], A, B)
    s <- score_groups(r$fc, r$ok, species_from_names(x$Protein.Names), HYE_TRUTH, centre = "Human")
    data.frame(set = label, species = s$group, n = s$n, median = s$median, IQR = s$IQR,
               within_0.5 = s$within, within_0.5_centred = s$within_centred)
  }, files, labels))
  bundled <- utils::read.csv(file.path(root, dir, ref), check.names = FALSE)
  rownames(res) <- NULL
  check(paste("re-creates", file.path(dir, ref)), isTRUE(all.equal(res, bundled, check.attributes = FALSE)))
}
recreate("diann-hye", c("report.pg_matrix.tsv", "pg_mbr.tsv", "pg_matrix_nombr_usequant.tsv"),
         c("pg_matrix (MBR, DIA-NN default filter)", "MBR pass, Q<=0.01 & PG.Q<=0.01", "no MBR (first-pass .quant reused)"),
         "diann_scores.csv")
recreate("diann-2.7-hye", c("report.pg_matrix.tsv", "pg_mbr.tsv"),
         c("pg_matrix (MBR, DIA-NN default filter)", "MBR pass, Q<=0.01 & PG.Q<=0.01"), "diann_scores_270.csv")

if (fails > 0) { cat(fails, "test(s) failed\n"); quit(status = 1) }
cat("all R tests passed\n")
