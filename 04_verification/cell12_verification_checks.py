# Cell 12 (add-on, run AFTER Run all in the same session): VERIFICATION CHECKS. Aggregates only.
say("\n=== CELL 12: verification checks ===")
# V1-V3: does the handling of the one admission with no baseline creatinine matter?
say("Admissions with no baseline creatinine:", n_((d.base24.isna()).sum()))
def rr_only(data, label):
    m = fit(data, 'aki', T(), label, 'poisson', quiet=True)
    if m is not None:
        r = RESULTS[-1]; say(f"{label}: N={r['N']}, RR {r['est']:.2f} ({r['lo']:.2f} to {r['hi']:.2f})")
rr_only(d[d.base24.notna()], "V1 Primary model excluding it")
for val in (0, 1):
    dv = d.copy(); dv.loc[dv.base24.isna(), 'aki'] = val
    rr_only(dv, f"V{2+val} Primary model with it counted as AKI={val}")
say("(If V1, V2 and V3 all round to 1.87 or differ only in the second decimal, the handling of that admission is immaterial.)")

# V4: independent recount of 'no creatinine before the procedure day' (different code path from Cell 3)
L = LABS[(LABS.itemid == CR_ITEM) & LABS.hadm_id.isin(ids)].dropna(subset=['valuenum']).merge(ADM[['hadm_id', 'admittime', 'dischtime']], on='hadm_id')
L = L[L.charttime <= L.dischtime]
w24 = L[L.charttime <= L.admittime + pd.Timedelta(hours=24)].merge(proc_date.reset_index(), on='hadm_id', how='left')
has_base = set(w24.hadm_id)
has_pre = set(w24[w24.charttime < w24.proc_date].hadm_id)
indep = len(has_base - has_pre)
orig = int((d.base24.notna() & d.base_pre.isna()).sum())
say(f"V4 Admissions with a baseline but none before the procedure day: independent recount {indep}, Cell 3 {orig}, match: {indep == orig}; admissions with a baseline: {len(has_base)}")

# V5: independent rebuild of the R3 (ICU-arrival) subset
if icu_path:
    i2 = pd.read_csv(icu_path, usecols=['hadm_id', 'intime'], parse_dates=['intime'])
    i2 = i2[i2.hadm_id.isin(ids)].groupby('hadm_id').intime.min().rename('icu_first').reset_index()
    X = L.merge(i2, on='hadm_id')
    cut = X[['icu_first']].assign(a=X.admittime + pd.Timedelta(hours=24)).min(axis=1)
    b_ok = set(X[X.charttime <= cut].hadm_id); p_ok = set(X[X.charttime > X.icu_first].hadm_id)
    r3 = b_ok & p_ok
    say(f"V5 R3 subset: independent rebuild {len(r3)}, Cell 3 {len(d2)}, same admissions: {r3 == set(d2.hadm_id)}")
else:
    say("V5 skipped: icustays file not found.")
