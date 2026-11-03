"""Shared scoring for known-answer benchmarks (numpy only). Python twin of common/known_ratio.R.

    import sys; sys.path.insert(0, "<repo>/common"); import known_ratio as kr

The rules match the existing HYE scripts: common/tests/test_known_ratio.py re-creates diann-hye/diann_scores.csv
with these functions, exactly as the R tests do.
"""
import csv
import re

import numpy as np

# PXD028735 (HYE): true log2(B/A) per species, and how to recognise each species
HYE_TRUTH = {"Human": 0.0, "Yeast": -1.0, "Ecoli": 2.0}
HYE_NAME_SUFFIX = {"Human": "_HUMAN", "Yeast": "_YEAST", "Ecoli": "_ECOLI"}   # UniProt entry names (DIA-NN, FragPipe)
HYE_TAXON = {"Human": "9606", "Yeast": "559292", "Ecoli": "83333"}            # OX= in FASTA headers (MaxQuant)

# Columns of the standard scores.csv (CONTRIBUTING.md)
SCORES_COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]


def species_from_names(names, suffixes=HYE_NAME_SUFFIX):
    """One label per name: the single group whose suffix occurs as a whole word; None if none or several."""
    pats = {g: re.compile(re.escape(p) + r"\b") for g, p in suffixes.items()}
    out = []
    for s in names:
        hits = [g for g, p in pats.items() if s is not None and p.search(str(s))]
        out.append(hits[0] if len(hits) == 1 else None)
    return out


def species_from_taxa(headers, taxa=HYE_TAXON):
    """Same from taxon IDs ("OX=9606"); IDs not in `taxa` are ignored."""
    by_id = {v: k for k, v in taxa.items()}
    out = []
    for h in headers:
        s = {by_id[i] for i in re.findall(r"OX=([0-9]+)", h or "") if i in by_id}
        out.append(s.pop() if len(s) == 1 else None)
    return out


def _log2(m):
    m = np.array(m, dtype=float)
    m[~np.isfinite(m) | (m <= 0)] = np.nan
    return np.log2(m)


def _row_means(x):
    """Mean of the non-missing values per row, summed left to right like R's rowMeans (NaN if none)."""
    valid = ~np.isnan(x)
    total = np.zeros(x.shape[0])
    for j in range(x.shape[1]):
        total = total + np.where(valid[:, j], x[:, j], 0.0)
    cnt = valid.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, total / np.maximum(cnt, 1), np.nan)


def log2_ratio(m, a, b, min_valid=2):
    """log2 ratio b / a per row from linear values (2-D array; a, b = column indices).
    Zeros, negatives and non-finite values are missing. Returns (ok, fc): ok = at least min_valid values in each
    group; fc = mean log2(b) - mean log2(a)."""
    lg = _log2(m)
    la, lb = lg[:, list(a)], lg[:, list(b)]
    ok = ((~np.isnan(la)).sum(axis=1) >= min_valid) & ((~np.isnan(lb)).sum(axis=1) >= min_valid)
    return ok, _row_means(lb) - _row_means(la)


def group_log2_means(m, groups, min_valid=2):
    """{group: mean log2 per row} for a dict of column-index lists; NaN where fewer than min_valid values."""
    lg = _log2(m)
    out = {}
    for g, cols in groups.items():
        x = lg[:, list(cols)]
        mu = _row_means(x)
        mu[(~np.isnan(x)).sum(axis=1) < min_valid] = np.nan
        out[g] = mu
    return out


def _r_round(x, digits):
    return float(np.round(x, digits)) if np.isfinite(x) else float("nan")


def score_groups(fc, ok, group, truth, tol=0.5, centre=None, digits=3):
    """Per group with a known truth: n, median and IQR of fc, share within tol of the truth (uncentred) and,
    with centre (e.g. "Human"), after subtracting that group's median. Returns (rows, centre_median)."""
    fc, ok, group = np.asarray(fc, dtype=float), np.asarray(ok, dtype=bool), np.asarray(group, dtype=object)
    shift = float(np.median(fc[ok & (group == centre)])) if centre is not None else float("nan")
    rows = []
    for g, e in truth.items():
        f = fc[ok & (group == g)]
        r = {"group": g, "n": int(f.size)}
        if f.size:
            q1, q3 = np.percentile(f, [25, 75])                     # same as R's default quantile type 7
            r.update(median=_r_round(np.median(f), digits), IQR=_r_round(q3 - q1, digits),
                     within=_r_round(np.mean(np.abs(f - e) <= tol), digits))
            if centre is not None:
                r["within_centred"] = _r_round(np.mean(np.abs(f - shift - e) <= tol), digits)
        else:
            r.update(median=float("nan"), IQR=float("nan"), within=float("nan"))
            if centre is not None:
                r["within_centred"] = float("nan")
        rows.append(r)
    return rows, shift


def call_outcome(called, fc, group, truth):
    """Per call: "right", "wrong_direction", "false" (truth 0), "not_called" or "no_truth"."""
    out = []
    for c, f, g in zip(called, fc, group):
        if not c:
            out.append("not_called")
        elif g not in truth:
            out.append("no_truth")
        elif truth[g] == 0:
            out.append("false")
        else:
            out.append("right" if np.sign(f) == np.sign(truth[g]) else "wrong_direction")
    return out


def call_summary(outcome):
    """Counts of call outcomes; fdp = (false + wrong_direction) / calls with a truth."""
    n = {k: sum(1 for o in outcome if o == k) for k in ("right", "wrong_direction", "false", "no_truth")}
    judged = n["right"] + n["wrong_direction"] + n["false"]
    return {"calls": judged + n["no_truth"], "right": n["right"], "wrong_direction": n["wrong_direction"],
            "false": n["false"], "fdp": (n["false"] + n["wrong_direction"]) / max(judged, 1)}


def write_scores(rows, path):
    """rows: list of dicts with exactly the SCORES_COLUMNS keys."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SCORES_COLUMNS)
        w.writeheader()
        for r in rows:
            if set(r) != set(SCORES_COLUMNS):
                raise ValueError(f"scores row has keys {sorted(r)}, expected {SCORES_COLUMNS}")
            w.writerow(r)
