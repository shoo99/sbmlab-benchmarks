import numpy as np, pandas as pd
rng=np.random.default_rng(0)
protein_matrix=pd.DataFrame(rng.normal(size=(60,40)))
labels=np.r_[np.zeros(30),np.ones(30)]
protein_matrix.iloc[30:, :5]+=1.0
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import numpy as np

# Load proteomics data
X = protein_matrix.values
y = labels

# Nested cross-validation: the inner loop tunes, the outer loop estimates performance
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

outer_scores = []
for train_idx, test_idx in outer_cv.split(X, y):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # Scaling lives inside the pipeline, so it is refit within every fold (no leakage)
    pipe = make_pipeline(StandardScaler(), RandomForestClassifier(random_state=42))
    search = GridSearchCV(pipe, {"randomforestclassifier__max_features": ["sqrt", 0.1]}, cv=inner_cv)
    search.fit(X_train, y_train)

    outer_scores.append(search.score(X_test, y_test))

print(f"CV accuracy: {np.mean(outer_scores):.3f} ± {np.std(outer_scores):.3f}")

pd.DataFrame({'Q.Value':rng.uniform(0,0.02,500),'Protein.Group':rng.choice(['P1','P2','P3'],500),'Precursor.Quantity':rng.lognormal(10,1,500)}).to_csv('huge_proteomics_dataset.csv',index=False)
import pandas as pd
import polars as pl  # multithreaded, lazy; see the timing below

# For very large proteomics tables (millions of rows)

df = pl.scan_csv("huge_proteomics_dataset.csv")
result = (
    df
    .filter(pl.col("Q.Value") < 0.01)
    .group_by("Protein.Group")
    .agg([
        pl.col("Precursor.Quantity").mean().alias("mean_intensity"),
        pl.col("Precursor.Quantity").std().alias("std_intensity"),
        pl.len().alias("n_peptides")
    ])
    .collect()
)

print(result)
