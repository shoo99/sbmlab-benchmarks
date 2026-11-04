"""GO over-representation helpers for go-release-drift (issue #19): OBO and GAF readers, propagation, the
hypergeometric upper tail and Benjamini-Hochberg. Standard library only; tested in test/test_go_ora.py."""
import gzip
import math


def read_obo(path):
    """Terms of one go-basic.obo: name, namespace, parents (is_a, part_of), alt ids, obsolete flag."""
    terms, alt, cur = {}, {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("["):
            cur = None
            if line == "[Term]":
                cur = {"id": None, "name": "", "ns": "", "parents": set(), "obsolete": False, "alt": []}
            continue
        if cur is None or ": " not in line:
            if line == "" and cur is not None and cur["id"]:
                terms[cur["id"]] = cur
            continue
        key, val = line.split(": ", 1)
        if key == "id":
            cur["id"] = val
        elif key == "name":
            cur["name"] = val
        elif key == "namespace":
            cur["ns"] = val
        elif key == "is_a":
            cur["parents"].add(val.split()[0])
        elif key == "relationship" and val.startswith("part_of "):
            cur["parents"].add(val.split()[1])
        elif key == "alt_id":
            cur["alt"].append(val)
        elif key == "is_obsolete" and val == "true":
            cur["obsolete"] = True
    if cur is not None and cur["id"]:
        terms[cur["id"]] = cur
    for t in terms.values():
        for a in t["alt"]:
            alt[a] = t["id"]
    return terms, alt


def ancestors_fn(terms):
    memo = {}

    def anc(t):
        if t not in memo:
            out = {t}
            for p in terms[t]["parents"]:
                if p in terms:
                    out |= anc(p)
            memo[t] = frozenset(out)
        return memo[t]
    return anc


def read_gaf(path, terms, alt, evidence):
    """Direct annotations per UniProt accession for each evidence setting (a name -> keep(evidence code) dict).
    Keeps human UniProtKB rows without a NOT qualifier; maps alternative ids to the primary id; drops terms that are
    obsolete or missing in the ontology; strips isoform suffixes from accessions."""
    direct = {e: {} for e in evidence}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or f[0] != "UniProtKB" or "taxon:9606" not in f[12]:
                continue
            if "NOT" in f[3].split("|"):
                continue
            go = alt.get(f[4], f[4])
            if go not in terms or terms[go]["obsolete"]:
                continue
            acc = f[1].split("-")[0]
            for e, keep in evidence.items():
                if keep(f[6]):
                    direct[e].setdefault(acc, set()).add(go)
    return direct


def hyper_upper(k, n, K, N, lf):
    """P(X >= k) for X ~ hypergeometric(N, K successes, n draws); lf = log-factorial table."""
    if k <= 0:
        return 1.0
    hi = min(K, n)
    if k > hi:
        return 0.0

    def logc(a, b):
        return lf[a] - lf[b] - lf[a - b]
    logp = logc(K, k) + logc(N - K, n - k) - logc(N, n)
    term, total, i = 1.0, 1.0, k
    mode = (n + 1) * (K + 1) // (N + 2)
    while i < hi:
        term *= (K - i) * (n - i) / ((i + 1) * (N - K - n + i + 1))
        total += term
        i += 1
        if i > mode and term < total * 1e-16:
            break
    return min(1.0, math.exp(logp) * total)


def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q, prev = [1.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def log_factorials(n):
    lf = [0.0]
    for i in range(1, n + 1):
        lf.append(lf[-1] + math.log(i))
    return lf
