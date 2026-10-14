"""Bader et al. 2020 (Mol Syst Biol, PMC7266499) Dataset EV3 (log10 protein intensities, 1,484 protein groups x 197 CSF samples)
+ EV2 (sample info). Writes a tidy matrix and labels: AD (biochemical) = 1, control = 0."""
import sys, pandas as pd
d = sys.argv[1]
e = pd.read_excel(f"{d}/MSB-16-e9356-s005.xlsx", sheet_name=0, header=None)
samples = e.iloc[1, 3:].tolist()
X = e.iloc[2:, 3:].astype(float); X.columns = samples
X.index = e.iloc[2:, 1].fillna(e.iloc[2:, 0]).astype(str)          # gene name, else accession
m = pd.read_excel(f"{d}/MSB-16-e9356-s004.xlsx", sheet_name=0)
m = m[m["primary biochemical AD classification"].notna()].set_index("sample name")
m = m.loc[[s for s in samples if s in m.index]]
X = X[m.index].T
y = (m["primary biochemical AD classification"] == "biochemical AD").astype(int)
X.to_csv("X_log10.csv"); pd.DataFrame({"y": y, "cohort": m["collection site"]}).to_csv("labels.csv")
print(X.shape, y.value_counts().to_dict(), m["collection site"].value_counts().to_dict())
print("missing cells: %.1f%%" % (100 * X.isna().mean().mean()))
