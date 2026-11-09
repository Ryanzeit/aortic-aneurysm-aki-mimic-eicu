# Cell 16 (add-on, run AFTER Run all in the same session). CREATININE RANGE AND DUPLICATE CHECK. DESCRIPTIVE ONLY: no outcome model, no new estimate.
# Why: RECORD 12.2 (data cleaning). The manuscript says creatinine was used as recorded, without plausibility limits or duplicate removal.
# This cell only looks at the values. Aggregates only; counts under SMALL_CELL print as '<11'.
say("\n=== CELL 16: creatinine range and duplicate check (descriptive; no outcome analysis) ===")
v = cr.valuenum
say("Creatinine values used in the analysis (item", CR_ITEM, "charted up to discharge, cohort admissions):", n_(len(v)), "| admissions with at least one:", n_(cr.hadm_id.nunique()))
say("Minimum / 1st / 50th / 99th / 99.9th percentile / maximum, mg/dL:", v.min(), "/", *[round(float(x), 2) for x in v.quantile([.01, .5, .99, .999])], "/", v.max())
for thr in [0.1, 0.2]:
    say(f"Values at or below {thr} mg/dL:", n_((v <= thr).sum()))
for thr in [10, 15, 20, 30]:
    say(f"Values at or above {thr} mg/dL:", n_((v >= thr).sum()), "| in how many admissions:", n_(cr[v >= thr].hadm_id.nunique()))
say("Zero or negative values:", n_((v <= 0).sum()), "| missing numeric value (should be 0 after cohort build):", n_(cr.valuenum.isna().sum()))
dup = cr.duplicated(subset=['hadm_id', 'charttime', 'valuenum'], keep=False)
say("Rows that exactly duplicate another row (same admission, charttime and value):", n_(dup.sum()))
same_time = cr.duplicated(subset=['hadm_id', 'charttime'], keep=False) & ~dup
say("Rows sharing an admission and charttime with a DIFFERENT value:", n_(same_time.sum()))
pk = d.peak_all.dropna(); bs = d.base24.dropna()
say("Admission-level peak creatinine, mg/dL (min / median / max):", round(float(pk.min()), 2), "/", round(float(pk.median()), 2), "/", round(float(pk.max()), 2),
    "| baseline (first 24 h) min / median / max:", round(float(bs.min()), 2), "/", round(float(bs.median()), 2), "/", round(float(bs.max()), 2))
say("Admissions with peak >= 15 mg/dL:", n_((pk >= 15).sum()), "| baseline <= 0.3 mg/dL:", n_((bs <= 0.3).sum()))
