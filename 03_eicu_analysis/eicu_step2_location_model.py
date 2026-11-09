# eICU step 2: location (thoracic vs abdominal) and AKI in eICU-CRD, with a text-derived comorbidity count (diagnosis table UNION pastHistory).
# Patterns below were fixed from the step-1 vocabulary BEFORE any outcome was examined. Aggregates only; counts under 11 are masked.
# eICU cannot identify open vs endovascular repair, so this answers the location question only.
import warnings, re, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
try:
    from google.colab import drive; drive.mount('/content/drive')
except Exception: pass
P = "/content/drive/MyDrive/eICU-CRD/"; SMALL = 11
def n_(v):
    v = int(v); return f"<{SMALL}" if 0 < v < SMALL else str(v)

allc = pd.read_csv(P+"aortic_aneurysm_eicu_cohort_final_v2.csv")
fin = set(pd.read_csv(P+"eicu_final_cohort_ids.csv").iloc[:, 0])
coh = allc[allc.patientunitstayid.isin(fin)].copy(); ids = set(coh.patientunitstayid)
print(f"cohort file rows {len(allc)} (events {n_(allc.aki.sum())}) | final-ID file stays {len(fin)} | analysed stays {len(coh)} | AKI events {n_(coh.aki.sum())}  (expect 382 stays, 86 events)")
print("final IDs missing from the cohort file:", n_(len(fin - set(allc.patientunitstayid))), "| location counts:", {k: n_(v) for k, v in coh.location.value_counts().items()})
print("missing values in model columns:", {c: int(coh[c].isna().sum()) for c in ['aki','location','age_numeric','gender','rupture','vasopressor','cci']})

# ---- text-derived comorbidity domains (seven Charlson-type domains + hypertension reported separately) ----
PH = {'htn': r'hypertension requiring treatment', 'dm': r'diabetes', 'ckd': r'renal insufficiency|chronic renal',
      'chf': r'congestive heart failure', 'copd': r'/copd/', 'cad': r'myocardial infarction|coronary artery bypass|procedural coronary intervention|coronary artery disease',
      'pvd': r'peripheral vascular disease', 'cvd': r'/strokes/|/tia\(s\)/'}
DX = {'htn': r'^hypertension', 'dm': r'diabetes mellitus', 'ckd': r'chronic renal|chronic kidney|renal insufficiency',
      'chf': r'^chf|congestive heart failure|heart failure', 'copd': r'^copd', 'cad': r'coronary artery disease|myocardial infarction|s/p cabg|s/p coronary',
      'pvd': r'peripheral vascular disease', 'cvd': r'stroke|\btia\b|cerebrovascular'}
ph = pd.read_csv(P+"pastHistory.csv.gz", usecols=['patientunitstayid','pasthistorypath']); ph = ph[ph.patientunitstayid.isin(ids)]
dx = pd.read_csv(P+"diagnosis.csv.gz", usecols=['patientunitstayid','diagnosisstring','icd9code']); dx = dx[dx.patientunitstayid.isin(ids)]
phl = ph.pasthistorypath.str.lower(); dxl = dx.diagnosisstring.str.split('|').str[-1].str.lower()
rows = []
for k in PH:
    a = set(ph[phl.str.contains(PH[k], regex=True, na=False)].patientunitstayid); b = set(dx[dxl.str.contains(DX[k], regex=True, na=False)].patientunitstayid)
    coh[k] = coh.patientunitstayid.isin(a | b).astype(int)
    rows.append((k, n_(len(a & ids)), n_(len(b & ids)), n_(len(a | b))))
print("\nComorbidity flags among analysed stays (distinct stays; small counts masked):")
print(pd.DataFrame(rows, columns=['domain','pastHistory','diagnosis','either']).to_string(index=False))
dom = ['dm','ckd','chf','copd','cad','pvd','cvd']; coh['comorb'] = coh[dom].sum(axis=1)
print("Text-derived comorbidity count (7 domains), distribution:", {int(k): n_(v) for k, v in coh.comorb.value_counts().sort_index().items()})
hasicd = set(dx[dx.icd9code.notna() & (dx.icd9code.astype(str).str.strip() != '')].patientunitstayid)
print(f"Stays with at least one diagnosis row carrying an ICD-9 code: {len(hasicd & ids)} of {len(coh)} | has_diagnosis_data==1 in the cohort file: {int((coh.has_diagnosis_data==1).sum())}")
try:
    pt = pd.read_csv(P+"patient.csv.gz", usecols=['patientunitstayid','hospitalid']); coh = coh.merge(pt, on='patientunitstayid', how='left')
except Exception as e: print("patient.csv.gz hospitalid not available:", type(e).__name__)

# ---- models: modified Poisson, robust SE ----
coh['thoracic'] = (coh.location == 'thoracic').astype(int)
def fit(data, f, label, cluster=None):
    try:
        kw = dict(cov_type='cluster', cov_kwds={'groups': data[cluster].values}) if cluster else dict(cov_type='HC1')
        m = smf.glm(f, data=data, family=sm.families.Poisson()).fit(**kw)
        b, se = m.params['thoracic'], m.bse['thoracic']; lo, hi = np.exp(b - 1.96*se), np.exp(b + 1.96*se)
        print(f"{label}: N={int(m.nobs)}, events={n_(data.aki.sum())}, thoracic vs abdominal RR {np.exp(b):.2f} ({lo:.2f} to {hi:.2f}), p={m.pvalues['thoracic']:.3f} | log RR {b:.4f}, SE {se:.4f}")
    except Exception as e: print(f"{label}: MODEL FAILED ({type(e).__name__}: {str(e)[:70]})")
print("\n=== eICU location models (post hoc) ===")
base = "aki ~ thoracic + age_numeric + C(gender) + rupture"
fit(coh, "aki ~ thoracic", "E0 unadjusted")
fit(coh, base + " + comorb", "E1 PRIMARY eICU: age, sex, rupture, text-derived comorbidity count")
fit(coh, base + " + comorb + htn", "E1b + hypertension")
fit(coh, base + " + cci", "E2 original diagnosis-table cci instead of the union count (old specification)")
fit(coh, base + " + comorb + vasopressor", "E3 + vasopressor (post-procedure; mediator risk)")
fit(coh[coh.rupture == 0], "aki ~ thoracic + age_numeric + C(gender) + comorb", "E4 excluding ruptured stays")
dissect = set(dx[dxl.str.contains('dissection', na=False)].patientunitstayid)
fit(coh[~coh.patientunitstayid.isin(dissect)], base + " + comorb", "E5 excluding stays with a dissection diagnosis")
if 'hospitalid' in coh and coh.hospitalid.notna().all():
    print("Hospitals in analysed stays:", n_(coh.hospitalid.nunique()))
    fit(coh, base + " + comorb", "E6 same as E1, SE clustered by hospital", cluster='hospitalid')
print("\nNext: the MIMIC location model (Cell 8) and the formal comparison.")
