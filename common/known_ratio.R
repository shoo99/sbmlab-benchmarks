## Shared scoring for known-answer benchmarks (base R only). See common/README.md and CONTRIBUTING.md.
##   source(file.path(<repo>, "common", "known_ratio.R"))
## The rules match the existing HYE scripts: common/tests/test_known_ratio.R re-creates diann-hye/diann_scores.csv
## with these functions.

## PXD028735 (HYE): true log2(B/A) per species, and how to recognise each species
HYE_TRUTH <- c(Human = 0, Yeast = -1, Ecoli = 2)
HYE_NAME_SUFFIX <- c(Human = "_HUMAN", Yeast = "_YEAST", Ecoli = "_ECOLI")     # UniProt entry names (DIA-NN, FragPipe)
HYE_TAXON <- c(Human = "9606", Yeast = "559292", Ecoli = "83333")             # OX= in FASTA headers (MaxQuant)

## Columns of the standard scores.csv (CONTRIBUTING.md)
SCORES_COLUMNS <- c("folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n")

## One group label per element: the single group whose name suffix occurs as a whole word; NA if none or several.
species_from_names <- function(x, suffixes = HYE_NAME_SUFFIX) {
  hits <- vapply(suffixes, function(p) grepl(paste0(p, "\\b"), x, perl = TRUE), logical(length(x)))
  hits <- matrix(hits, nrow = length(x))
  out <- rep(NA_character_, length(x))
  one <- rowSums(hits) == 1
  out[one] <- names(suffixes)[max.col(hits[one, , drop = FALSE], ties.method = "first")]
  out
}

## Same from taxon IDs ("OX=9606"); IDs not in `taxa` are ignored, as in maxquant-hye/score_quant.R.
species_from_taxa <- function(x, taxa = HYE_TAXON) {
  vapply(x, function(h) {
    ox <- unique(sub("^OX=", "", regmatches(h, gregexpr("OX=[0-9]+", h))[[1]]))
    s <- unique(names(taxa)[match(ox, taxa)])
    s <- s[!is.na(s)]
    if (length(s) == 1) s else NA_character_
  }, "", USE.NAMES = FALSE)
}

## log2 ratio b / a per row from linear values. Zeros, negatives and non-finite values are missing.
## ok: at least `min_valid` values in each group; fc: mean log2(b) - mean log2(a).
log2_ratio <- function(m, a, b, min_valid = 2) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  m[!is.finite(m) | m <= 0] <- NA
  l <- log2(m)
  ok <- rowSums(!is.na(l[, a, drop = FALSE])) >= min_valid & rowSums(!is.na(l[, b, drop = FALSE])) >= min_valid
  fc <- rowMeans(l[, b, drop = FALSE], na.rm = TRUE) - rowMeans(l[, a, drop = FALSE], na.rm = TRUE)
  list(ok = unname(ok), fc = unname(fc))
}

## Mean log2 value per group of columns (named list of column sets); NA where a group has fewer than min_valid values.
group_log2_means <- function(m, groups, min_valid = 2) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  m[!is.finite(m) | m <= 0] <- NA
  l <- log2(m)
  sapply(groups, function(cols) {
    x <- l[, cols, drop = FALSE]
    mu <- rowMeans(x, na.rm = TRUE)
    mu[rowSums(!is.na(x)) < min_valid] <- NA
    unname(mu)
  })
}

## Per group with a known truth: n, median and IQR of fc, and the share within `tol` of the truth (uncentred).
## centre: name of the group whose median is subtracted for within_centred (e.g. "Human"), or NULL.
score_groups <- function(fc, ok, group, truth, tol = 0.5, centre = NULL, digits = 3) {
  shift <- if (is.null(centre)) NA_real_ else median(fc[ok & group %in% centre])
  rows <- lapply(names(truth), function(g) {
    i <- ok & group %in% g
    f <- fc[i]
    e <- truth[[g]]
    r <- data.frame(group = g, n = sum(i), median = round(median(f), digits), IQR = round(IQR(f), digits),
                    within = round(mean(abs(f - e) <= tol), digits), stringsAsFactors = FALSE)
    if (!is.null(centre)) r$within_centred <- round(mean(abs(f - shift - e) <= tol), digits)
    r
  })
  out <- do.call(rbind, rows)
  attr(out, "centre_median") <- shift
  out
}

## Outcome of each call against the truth: "right" (changed and in the true direction), "wrong_direction",
## "false" (called although the truth is 0), "not_called", or "no_truth" (called, group without a truth).
call_outcome <- function(called, fc, group, truth) {
  e <- unname(truth[group])
  called <- !is.na(called) & called
  ifelse(!called, "not_called",
         ifelse(is.na(e), "no_truth",
                ifelse(e == 0, "false", ifelse(sign(fc) == sign(e), "right", "wrong_direction"))))
}

## Counts of call outcomes; fdp = (false + wrong_direction) / calls with a truth.
call_summary <- function(outcome) {
  n <- function(k) sum(outcome == k)
  judged <- n("right") + n("wrong_direction") + n("false")
  data.frame(calls = judged + n("no_truth"), right = n("right"), wrong_direction = n("wrong_direction"),
             false = n("false"), fdp = (n("false") + n("wrong_direction")) / max(judged, 1))
}

## Rows in the standard scores.csv format; arguments are recycled like data.frame().
scores_rows <- function(folder, dataset, tool, version, setting, subset, metric, value, n) {
  data.frame(folder = folder, dataset = dataset, tool = tool, version = version, setting = setting,
             subset = subset, metric = metric, value = value, n = n, stringsAsFactors = FALSE)
}

write_scores <- function(rows, path) {
  stopifnot(identical(names(rows), SCORES_COLUMNS))
  utils::write.csv(rows, path, row.names = FALSE)
}
