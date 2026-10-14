"""Train on Berlin + Sweden, test on Magdeburg + Kiel (never touched during training)."""
import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

X = pd.read_csv("X_log10.csv", index_col=0); lab = pd.read_csv("labels.csv", index_col=0)
y = lab["y"].values; train = lab["cohort"].isin(["Berlin", "Sweden"]).values; test = ~train
K = 20
print("train", train.sum(), "AD", y[train].sum(), "| test", test.sum(), "AD", y[test].sum())
pipe = lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                             SelectKBest(f_classif, k=K), LogisticRegression(max_iter=5000))

def external(y, leaky):
    if leaky:   # proteins chosen using every sample, test cohorts included
        Xi = SimpleImputer(strategy="median").fit_transform(X)
        cols = SelectKBest(f_classif, k=K).fit(Xi, y).get_support()
        m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=5000))
        m.fit(X.loc[train, cols], y[train]); return roc_auc_score(y[test], m.predict_proba(X.loc[test, cols])[:, 1])
    m = pipe().fit(X[train], y[train]); return roc_auc_score(y[test], m.predict_proba(X[test])[:, 1])

print("real labels, proper:", round(external(y, False), 3), "| leaky:", round(external(y, True), 3))
m = pipe().fit(X[train], y[train]); sel = X.columns[m.named_steps["selectkbest"].get_support()]
print("proteins selected on the training cohorts:", ", ".join(sel))
rng = np.random.default_rng(1); res = []
for i in range(50):
    ys = rng.permutation(y)
    res.append(dict(rep=i, proper=external(ys, False), leaky=external(ys, True)))
r = pd.DataFrame(res); r.to_csv("external_shuffled.csv", index=False)
print("shuffled labels, 50 permutations:"); print(r[["proper", "leaky"]].describe().loc[["mean", "std", "min", "max"]].round(3))
