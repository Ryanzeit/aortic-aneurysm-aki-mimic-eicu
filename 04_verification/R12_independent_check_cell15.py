# R12: INDEPENDENT RECOMPUTATION of Cell 12 AND the abstract-bound Cell 15 results. (It repeats the R11 checks first.) ORIGINAL R11 header follows:
# R11: INDEPENDENT RECOMPUTATION of Cell 12 (pre-admission baseline lookback). Colab cell. Aggregates only; counts below 11 are masked.
# Different code path from the notebook on purpose: DuckDB SQL for the lookback window, median, peak and AKI flags,
# and the same hand-written Poisson IRLS with a patient-clustered sandwich as R10 (numpy only; no statsmodels).
# What it checks: the 7-365 day window, the per-admission median, coverage by approach, the in-hospital/pre-admission ratio table,
#                 and P0, P1, P2 (point estimates and CIs).
# What it does NOT check: the 521 "no creatinine before the procedure day" figure (needs procedure dates), the exclusions behind the 740,
#                 or anything about other analyses.
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
LE  = P+"mimic-iv/hosp/labevents.csv.gz"   # large; reading takes several minutes. Edit if Cell 12 needed a different path.
CR_ITEM = 50912                              # serum creatinine; must match CR_ITEM in the notebook
LO_DAYS, HI_DAYS = 7, 365                    # the window Cell 12 uses; do NOT change it here
SMALL = 11
def m(v):
    v = int(v); return f"<{SMALL}" if 0 < v < SMALL else str(v)

con = duckdb.connect()
con.execute(f"CREATE VIEW coh AS SELECT * FROM read_csv_auto('{COH}', header=true)")
con.execute(f"""CREATE VIEW adm AS SELECT subject_id AS sid, hadm_id, CAST(admittime AS TIMESTAMP) AS a, CAST(dischtime AS TIMESTAMP) AS dsc,
                admission_type FROM read_csv_auto('{ADM}', header=true)""")
con.execute(f"""CREATE VIEW cr AS SELECT hadm_id, CAST(charttime AS TIMESTAMP) AS t, valuenum AS v FROM read_csv_auto('{CRF}', header=true)
                WHERE itemid = {CR_ITEM} AND valuenum IS NOT NULL""")
# creatinine for cohort patients from the full lab table, ALL dates (includes outpatient draws with no admission id)
con.execute(f"""CREATE TABLE le AS
  SELECT subject_id AS sid, CAST(charttime AS TIMESTAMP) AS t, valuenum AS v
  FROM read_csv('{LE}', header=true, types={{'itemid':'INTEGER','valuenum':'DOUBLE','charttime':'TIMESTAMP','subject_id':'BIGINT'}})
  WHERE itemid = {CR_ITEM} AND valuenum IS NOT NULL
    AND subject_id IN (SELECT DISTINCT a.sid FROM adm a JOIN coh c USING (hadm_id))""")
print("Creatinine values for cohort patients in labevents (all dates):", con.execute("SELECT count(*) FROM le").fetchone()[0],
      "| patients:", con.execute("SELECT count(DISTINCT sid) FROM le").fetchone()[0], "   (notebook: 24937 | 689)")

# per admission: first-24h baseline (as R10), peak up to discharge (as R10), pre-admission median 7-365 days before admission
con.execute(f"""CREATE TABLE d AS
 SELECT c.*, a.sid, a.admission_type, a.a AS admit_t,
   (SELECT min(v) FROM cr WHERE cr.hadm_id=c.hadm_id AND cr.t <= a.dsc AND cr.t <= a.a + INTERVAL 24 HOUR) AS base24,
   (SELECT max(v) FROM cr WHERE cr.hadm_id=c.hadm_id AND cr.t <= a.dsc) AS peak,
   (SELECT median(le.v) FROM le WHERE le.sid=a.sid AND le.t <= a.a - INTERVAL {LO_DAYS} DAY AND le.t >= a.a - INTERVAL {HI_DAYS} DAY) AS prior,
   (SELECT count(*) FROM le WHERE le.sid=a.sid AND le.t <= a.a - INTERVAL {LO_DAYS} DAY AND le.t >= a.a - INTERVAL {HI_DAYS} DAY) AS n_prior
 FROM coh c JOIN adm a USING (hadm_id)""")
con.execute("""CREATE TABLE a2 AS SELECT *,
   CASE WHEN prior IS NOT NULL THEN 1 ELSE 0 END AS has_prior,
   coalesce(prior, base24) AS base_hier,
   CASE WHEN prior IS NOT NULL AND peak/prior >= 1.5 THEN 1 ELSE 0 END AS aki_prior,
   CASE WHEN coalesce(prior, base24) IS NOT NULL AND peak/coalesce(prior, base24) >= 1.5 THEN 1 ELSE 0 END AS aki_hier
   FROM d""")
N = con.execute("SELECT count(*) FROM a2").fetchone()[0]
print(f"\nAdmissions: {N} (expect 740) | patients: {con.execute('SELECT count(DISTINCT sid) FROM a2').fetchone()[0]} (expect 689)")
hp = con.execute("SELECT sum(has_prior) FROM a2").fetchone()[0]
print(f"Admissions with a pre-admission baseline: {m(hp)} of {N} ({100*hp/N:.1f}%)   (notebook: 395 of 740)")
print("By approach (approach, has_prior, n); notebook: endovascular 236/243, open 109/152:")
for r in con.execute("SELECT approach, has_prior, count(*) FROM a2 GROUP BY 1,2 ORDER BY 1,2").fetchall(): print("  ", r[0], r[1], m(r[2]))
print("Prior values per admission with a baseline (median, IQR):",
      con.execute("SELECT median(n_prior), quantile_cont(n_prior,0.25), quantile_cont(n_prior,0.75) FROM a2 WHERE has_prior=1").fetchone())

# bias diagnostic: first-24h baseline divided by pre-admission baseline
print("\nIn-hospital baseline / pre-admission baseline by approach: 25%, 50%, 75% (notebook: endo 0.82/0.91/1.00, open 0.81/0.92/1.09)")
for r in con.execute("""SELECT approach, round(quantile_cont(base24/prior,0.25),2), round(quantile_cont(base24/prior,0.5),2), round(quantile_cont(base24/prior,0.75),2)
                        FROM a2 WHERE has_prior=1 AND base24 IS NOT NULL GROUP BY 1""").fetchall(): print("  ", r)
print("Counts with ratio >= 1.25 (hi) and <= 0.80 (lo); notebook: endo N243 hi21 lo58, open N152 hi19 lo38:")
for r in con.execute("""SELECT approach, count(*), sum(CASE WHEN base24/prior >= 1.25 THEN 1 ELSE 0 END), sum(CASE WHEN base24/prior <= 0.80 THEN 1 ELSE 0 END)
                        FROM a2 WHERE has_prior=1 AND base24 IS NOT NULL GROUP BY 1""").fetchall(): print("  ", r[0], "N", m(r[1]), "hi", m(r[2]), "lo", m(r[3]))

print("\nAKI events: manuscript flag (aki) total", m(con.execute("SELECT sum(aki) FROM a2").fetchone()[0]),
      "| in has_prior subset", m(con.execute("SELECT sum(aki) FROM a2 WHERE has_prior=1").fetchone()[0]), "(notebook P0: 110)")
print("aki_prior in subset", m(con.execute("SELECT sum(aki_prior) FROM a2").fetchone()[0]), "(notebook P1: 108) | aki_hier all", m(con.execute("SELECT sum(aki_hier) FROM a2").fetchone()[0]), "(notebook P2: 227)")
print("Admissions with no baseline at all (no prior and no 24 h value):", m(con.execute("SELECT count(*) FROM a2 WHERE base_hier IS NULL").fetchone()[0]))

df = con.execute("SELECT * FROM a2").df()
def urg(x):
    x = str(x).upper()
    if 'ELECTIVE' in x or 'SURGICAL SAME DAY' in x: return 'elective'
    if 'OBSERV' in x: return 'observation'
    return 'emergency_urgent'
df['urg'] = df.admission_type.map(urg)

def poisson_cluster(X, y, g):
    n, k = X.shape; b = np.zeros(k); b[0] = np.log(y.mean())
    for _ in range(100):
        mu = np.exp(X @ b)
        step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (y - mu)); b = b + step
        if np.abs(step).max() < 1e-10: break
    mu = np.exp(X @ b); bread = np.linalg.inv(X.T @ (X * mu[:, None]))
    s = X * (y - mu)[:, None]; ug, inv = np.unique(g, return_inverse=True)
    S = np.zeros((len(ug), k)); np.add.at(S, inv, s)
    G = len(ug); meat = S.T @ S * (G/(G-1)) * ((n-1)/(n-k))      # same small-sample factor statsmodels applies
    V = bread @ meat @ bread; return b, np.sqrt(np.diag(V)), mu
def design(sub, extra=(), unadj=False):
    def dum(col):                    # dummies for every level except the most common one (the reference does not change the open RR)
        vc = sub[col].astype(str).value_counts(); ref = vc.index[0]
        return [(f"{col}={l}", (sub[col].astype(str) == l).values.astype(float)) for l in sorted(vc.index) if l != ref]
    cols = [("const", np.ones(len(sub))), ("open", (sub.approach == 'open').values.astype(float))]
    if not unadj:
        cols += dum('location') + [("age", sub.age_at_admission.values.astype(float))] + dum('gender') + [("cci", sub.cci.values.astype(float))] \
              + dum('anchor_year_group') + [("rupture", sub.rupture.values.astype(float))] + dum('urg') \
              + [(e, sub[e].values.astype(float)) for e in extra]
    keep = [(n, c) for n, c in cols if n == "const" or c.std() > 0]          # drop columns with no variation in a subset (e.g. rupture when ruptured repairs are excluded)
    return np.column_stack([c for _, c in keep]), [n for n, _ in keep]
def logit_cluster(X, y, g):
    n, k = X.shape; b = np.zeros(k); b[0] = np.log(y.mean()/(1-y.mean()))
    for _ in range(100):
        mu = 1/(1+np.exp(-(X @ b))); W = mu*(1-mu)
        step = np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (y - mu)); b = b + step
        if np.abs(step).max() < 1e-10: break
    mu = 1/(1+np.exp(-(X @ b))); W = mu*(1-mu); bread = np.linalg.inv(X.T @ (X * W[:, None]))
    s = X * (y - mu)[:, None]; ug, inv = np.unique(g, return_inverse=True)
    S = np.zeros((len(ug), k)); np.add.at(S, inv, s)
    G = len(ug); meat = S.T @ S * (G/(G-1)) * ((n-1)/(n-k))
    V = bread @ meat @ bread; return b, np.sqrt(np.diag(V)), mu
from math import erfc, sqrt
def pfmt(z):
    p = erfc(abs(z)/sqrt(2)); return '<0.001' if p < 0.001 else f'{p:.3f}'
def run(label, sub, ycol, extra=(), unadj=False, logit=False):
    X, names = design(sub, extra, unadj); y = sub[ycol].values.astype(float)
    b, se, mu = (logit_cluster if logit else poisson_cluster)(X, y, sub.sid.values); j = names.index("open")
    print(f"{label}: N={len(y)}, events={m(int(y.sum()))}, {'OR' if logit else 'RR'} {np.exp(b[j]):.2f} ({np.exp(b[j]-1.96*se[j]):.2f} to {np.exp(b[j]+1.96*se[j]):.2f}), p={pfmt(b[j]/se[j])}")
    return X, names, b
df['ratio_h'] = df.peak / df.base_hier
df['stage2h'] = (df.ratio_h >= 2.0).astype(int); df['stage3h'] = (df.ratio_h >= 3.0).astype(int)
print("\nStage table on the hierarchy baseline (rows = approach; no AKI / stage 1 / stage 2 / stage 3); Cell 15: endovascular 358/48/29/43, open 154/45/37/25")
st = np.select([df.ratio_h >= 3.0, df.ratio_h >= 2.0, df.ratio_h >= 1.5, df.ratio_h < 1.5], [3, 2, 1, 0], default=-1)
for ap in ['endovascular', 'open']:
    print("  ", ap, [m(int(((st == k) & (df.approach == ap)).sum())) for k in [0, 1, 2, 3]], "| ratio missing or infinite:", m(int(((st == -1) & (df.approach == ap)).sum())))
print("\nCell 15 values to match:")
print("  H1 1.77 (1.41 to 2.22) 227 | H1b OR 2.71 (1.79 to 4.12) | H1u 1.64 (1.32 to 2.03) | H2 1.77 (1.27 to 2.47) 134 | H3 1.13 (0.67 to 1.90) 68")
print("  H4 1.84 (1.45 to 2.33) N689 | H5 1.76 (1.32 to 2.34) N623 | H6 1.56 (1.24 to 1.96) | crude 41.0% vs 25.1% | standardized 43.2% vs 24.1%, RD 19.1 | unadjusted RD 15.9")
run("R12 H1  primary, hierarchy baseline", df, 'aki_hier')
run("R12 H1b logistic odds ratio", df, 'aki_hier', logit=True)
run("R12 H1u unadjusted RR", df, 'aki_hier', unadj=True)
o = (df.approach.values == 'open'); yh = df.aki_hier.values
print(f"R12 crude AKI: open {100*yh[o].mean():.1f}% | endovascular {100*yh[~o].mean():.1f}% | unadjusted risk difference, point estimate {100*(yh[o].mean()-yh[~o].mean()):.1f} points")
X, names, b = run("R12 (logistic model for standardization)", df, 'aki_hier', logit=True)
j = names.index("open"); X1 = X.copy(); X1[:, j] = 1; X0 = X.copy(); X0[:, j] = 0
p1 = (1/(1+np.exp(-(X1 @ b)))).mean(); p0 = (1/(1+np.exp(-(X0 @ b)))).mean()
print(f"R12 standardized (g-computation, logistic): all open {100*p1:.1f}% vs all endovascular {100*p0:.1f}%; risk difference {100*(p1-p0):.1f} points; standardized RR {p1/p0:.2f}  (point estimates only; the bootstrap CI is not repeated)")
run("R12 H2  stage 2 or 3", df, 'stage2h')
run("R12 H3  stage 3", df, 'stage3h')
run("R12 H4  first admission per patient", df.sort_values('admit_t').drop_duplicates('sid'), 'aki_hier')
run("R12 H5  excluding ruptured repairs", df[df.rupture == 0], 'aki_hier')
run("R12 H6  add any-vasopressor", df, 'aki_hier', extra=('vasopressor',))
print("\nAgreement rule: counts must match exactly; RRs and ORs to 2 decimals; percentages to 0.1.")
print("NOT checked here: H7 (draws per day), H8 (overlap weights), H9/H10 (location strata and interaction), the bootstrap intervals, the E-value (formula only), the 298-of-521 figure, the 795 -> 740 exclusions.")
