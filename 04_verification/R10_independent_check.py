# R10: INDEPENDENT RECOMPUTATION of the primary result. Colab cell. Aggregates only; counts below 11 are masked.
# Different code path from the analysis notebook on purpose: DuckDB SQL for baseline/peak/AKI and cohort counts, and a hand-written
# Poisson IRLS with a patient-clustered sandwich (numpy only; no statsmodels, no pandas model formulas).
# What it can check: the creatinine rule, the AKI flag, covariate coding, and the primary and unadjusted estimates.
# What it cannot check: the exclusions that produced the 740 (that needs the original cohort-build tables; see the note at the end).
import subprocess, sys
try: import duckdb
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "duckdb"], check=True); import duckdb
import numpy as np
try:
    from google.colab import drive; drive.mount('/content/drive')
except Exception: pass

P = "/content/drive/MyDrive/"            # same as the notebook; edit only if the notebook needed edits
COH = P+"mimic-iv/aortic_aneurysm_final_modeling_dataset_v4.csv"
ADM = P+"mimic-iv/hosp/admissions.csv.gz"
CRF = P+"mimic-iv/aortic_aneurysm_creatinine_raw.csv"
SMALL = 11
def m(v):
    v = int(v); return f"<{SMALL}" if 0 < v < SMALL else str(v)

con = duckdb.connect()
con.execute(f"""CREATE VIEW coh AS SELECT * FROM read_csv_auto('{COH}', header=true)""")
con.execute(f"""CREATE VIEW adm AS SELECT subject_id AS sid, hadm_id, CAST(admittime AS TIMESTAMP) AS a, CAST(dischtime AS TIMESTAMP) AS dsc,
                admission_type FROM read_csv_auto('{ADM}', header=true)""")
con.execute(f"""CREATE VIEW cr AS SELECT hadm_id, CAST(charttime AS TIMESTAMP) AS t, valuenum AS v FROM read_csv_auto('{CRF}', header=true)
                WHERE itemid = 50912 AND valuenum IS NOT NULL""")
# baseline = lowest creatinine up to 24 h after admission, INCLUDING draws before the admission timestamp; peak = highest up to discharge
con.execute("""CREATE TABLE d AS
 SELECT c.*, a.sid, a.admission_type,
   (SELECT min(v) FROM cr WHERE cr.hadm_id=c.hadm_id AND cr.t <= a.dsc AND cr.t <= a.a + INTERVAL 24 HOUR) AS base,
   (SELECT max(v) FROM cr WHERE cr.hadm_id=c.hadm_id AND cr.t <= a.dsc) AS peak
 FROM coh c JOIN adm a USING (hadm_id)""")
n_rows = con.execute("SELECT count(*) FROM d").fetchone()[0]; n_coh = con.execute("SELECT count(*) FROM coh").fetchone()[0]
n_pat = con.execute("SELECT count(DISTINCT sid) FROM d").fetchone()[0]
print(f"Rows in cohort file {n_coh} | joined to admissions {n_rows} | patients {n_pat}")
print("approach counts:", {r[0]: m(r[1]) for r in con.execute("SELECT approach, count(*) FROM d GROUP BY 1").fetchall()})
nb = con.execute("SELECT count(*) FROM d WHERE base IS NOT NULL").fetchone()[0]
print(f"Admissions with a baseline: {nb} | without: {m(n_rows-nb)}")
con.execute("CREATE TABLE a2 AS SELECT *, CASE WHEN base IS NOT NULL AND peak/base >= 1.5 THEN 1 ELSE 0 END AS aki2 FROM d")
print("AKI events recomputed:", m(con.execute("SELECT sum(aki2) FROM a2").fetchone()[0]),
      "| manuscript flag:", m(con.execute("SELECT sum(aki) FROM a2").fetchone()[0]),
      "| disagreements:", m(con.execute("SELECT count(*) FROM a2 WHERE aki2 <> aki").fetchone()[0]))
print("Disagreements among admissions that HAVE a baseline:", m(con.execute("SELECT count(*) FROM a2 WHERE base IS NOT NULL AND aki2 <> aki").fetchone()[0]))
print("Admissions with NO baseline but flag aki=1:", m(con.execute("SELECT count(*) FROM a2 WHERE base IS NULL AND aki = 1").fetchone()[0]))

df = con.execute("SELECT * FROM a2").df()
def urg(x):
    x = str(x).upper()
    if 'ELECTIVE' in x or 'SURGICAL SAME DAY' in x: return 'elective'
    if 'OBSERV' in x: return 'observation'
    return 'emergency_urgent'
df['urg'] = df.admission_type.map(urg)
y = df['aki'].values.astype(float)            # the manuscript flag, as the notebook uses

def dummies(col, ref):
    lv = [l for l in sorted(df[col].astype(str).unique()) if l != ref]
    return [(f"{col}={l}", (df[col].astype(str) == l).values.astype(float)) for l in lv]
def design(terms):
    cols = [("const", np.ones(len(df)))] + terms
    return np.column_stack([c for _, c in cols]), [n for n, _ in cols]
open_ = [("open", (df.approach == 'open').values.astype(float))]
adj = (open_ + dummies('location', 'abdominal') + [("age", df.age_at_admission.values.astype(float))] + dummies('gender', sorted(df.gender.astype(str).unique())[0])
       + [("cci", df.cci.values.astype(float))] + dummies('anchor_year_group', sorted(df.anchor_year_group.astype(str).unique())[0])
       + [("rupture", df.rupture.values.astype(float))] + dummies('urg', 'elective'))

def poisson_cluster(X, y, g):
    n, k = X.shape; b = np.zeros(k); b[0] = np.log(y.mean())
    for _ in range(100):
        mu = np.exp(X @ b); W = mu
        step = np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (y - mu)); b = b + step
        if np.abs(step).max() < 1e-10: break
    mu = np.exp(X @ b); bread = np.linalg.inv(X.T @ (X * mu[:, None]))
    s = X * (y - mu)[:, None]; ug, inv = np.unique(g, return_inverse=True)
    S = np.zeros((len(ug), k)); np.add.at(S, inv, s)
    G = len(ug); meat = S.T @ S * (G/(G-1)) * ((n-1)/(n-k))      # same small-sample factor statsmodels applies
    V = bread @ meat @ bread; return b, np.sqrt(np.diag(V)), mu
def rr(label, terms):
    X, names = design(terms); b, se, mu = poisson_cluster(X, y, df.sid.values)
    j = names.index("open"); print(f"{label}: N={len(y)}, events={m(int(y.sum()))}, RR {np.exp(b[j]):.2f} ({np.exp(b[j]-1.96*se[j]):.2f} to {np.exp(b[j]+1.96*se[j]):.2f})")
    return X, b, mu
print("\nExpect from the notebook (real run): adjusted RR 1.87 (1.49-2.34); unadjusted RR 1.70 (1.37-2.10); RD 21.1 (adj) and 17.3 (unadj) points")
X, b, mu = rr("R10 adjusted (modified Poisson, cluster = patient)", adj)
rr("R10 unadjusted", open_)
# standardized risk difference (g-computation; point estimate only here)
j = 1; X1 = X.copy(); X1[:, j] = 1; X0 = X.copy(); X0[:, j] = 0
rd = (np.exp(X1 @ b).mean() - np.exp(X0 @ b).mean()) * 100
print(f"R10 standardized (adjusted) risk difference, point estimate: {rd:.1f} percentage points")
o=df.approach.values=='open'; print(f"R10 crude risk difference, point estimate: {100*(y[o].mean()-y[~o].mean()):.1f} percentage points")
print("\nNOTE: this does not re-derive the 795->740 exclusions. That needs icustays, procedures_icd, diagnoses_icd and the full creatinine table; a separate check would be needed for that.")
