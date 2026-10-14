"""Leaky vs proper biomarker pipelines on the Bader 2020 CSF proteome, with real and shuffled labels."""
import numpy as np, pandas as pd, sys
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score

X = pd.read_csv("X_log10.csv", index_col=0); lab = pd.read_csv("labels.csv", index_col=0)
y = lab["y"].values; K = int(sys.argv[1]) if len(sys.argv) > 1 else 20

def model():
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=5000))

def leaky_cv(X, y, seed):
    # selection on ALL samples (including those later used for testing), then CV only the classifier
    Xi = SimpleImputer(strategy="median").fit_transform(X)
    cols = SelectKBest(f_classif, k=K).fit(Xi, y).get_support()
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    return cross_val_score(model(), X.loc[:, cols], y, cv=cv, scoring="roc_auc").mean()

def proper_cv(X, y, seed):
    # selection inside the pipeline, so it is refitted on each training fold only
    pipe = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         SelectKBest(f_classif, k=K), LogisticRegression(max_iter=5000))
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    return cross_val_score(pipe, X, y, cv=cv, scoring="roc_auc").mean()

def train_auc(X, y):
    # no CV at all: select, fit and score on the same samples
    pipe = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         SelectKBest(f_classif, k=K), LogisticRegression(max_iter=5000))
    return roc_auc_score(y, pipe.fit(X, y).predict_proba(X)[:, 1])

rows = []
rng = np.random.default_rng(0)
for rep in range(50):
    ys = rng.permutation(y)                      # shuffled labels: true AUC = 0.5
    rows += [dict(labels="real", method=m, rep=rep, auc=f(X, y, rep)) for m, f in (("leaky CV", leaky_cv), ("proper CV", proper_cv))]
    rows += [dict(labels="shuffled", method=m, rep=rep, auc=f(X, ys, rep)) for m, f in (("leaky CV", leaky_cv), ("proper CV", proper_cv))]
    rows.append(dict(labels="shuffled", method="no CV (training AUC)", rep=rep, auc=train_auc(X, ys)))
rows.append(dict(labels="real", method="no CV (training AUC)", rep=0, auc=train_auc(X, y)))
r = pd.DataFrame(rows); r.to_csv(f"leak_k{K}.csv", index=False)
print(f"k = {K} proteins selected, 5-fold CV, 50 repeats (shuffled: 50 different permutations)")
print(r.groupby(["labels", "method"]).auc.describe()[["mean", "std", "min", "max"]].round(3))
