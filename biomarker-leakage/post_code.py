import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

# Bader et al. 2020, Dataset EV3 (log10 intensities) and EV2 (sample information)
e = pd.read_excel("MSB-16-e9356-s005.xlsx", header=None)
X = e.iloc[2:, 3:].astype(float).T
X.index = e.iloc[1, 3:]; X.columns = e.iloc[2:, 0]          # UniProt accessions (unique)
info = pd.read_excel("MSB-16-e9356-s004.xlsx").set_index("sample name")
info = info[info["primary biochemical AD classification"].notna()]
X = X.loc[[s for s in X.index if s in info.index]]
y = (info.loc[X.index, "primary biochemical AD classification"] == "biochemical AD").astype(int).values

cv = StratifiedKFold(5, shuffle=True, random_state=0)
clf = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=5000))
pipe = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                     SelectKBest(f_classif, k=20), LogisticRegression(max_iter=5000))
Xi = SimpleImputer(strategy="median").fit_transform(X)
rng = np.random.default_rng(0)

leaky, proper = [], []
for _ in range(20):
    y_null = rng.permutation(y)                          # shuffled labels: the honest answer is AUC 0.5

    # Leaky: pick the 20 "best" proteins using ALL samples, then cross-validate the classifier
    keep = SelectKBest(f_classif, k=20).fit(Xi, y_null).get_support()
    leaky.append(cross_val_score(clf, X.loc[:, keep], y_null, cv=cv, scoring="roc_auc").mean())

    # Proper: selection is a step of the pipeline, so it only ever sees the training folds
    proper.append(cross_val_score(pipe, X, y_null, cv=cv, scoring="roc_auc").mean())

print(X.shape, f"leaky AUC {np.mean(leaky):.2f} (max {max(leaky):.2f}) | proper AUC {np.mean(proper):.2f} (max {max(proper):.2f})")
