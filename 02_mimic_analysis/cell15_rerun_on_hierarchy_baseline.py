# Cell 15 (add-on, run AFTER Run all AND Cell 12 in the same session). MAIN-TEXT ANALYSES RE-RUN WITH THE HIERARCHY BASELINE. POST HOC, requested by Ryan 4 Nov 2026.
# Why: Cell 12 / R11 verified the hierarchy outcome (aki_hier: median creatinine 7-365 days before admission where it exists, otherwise the lowest value in the first 24 h).
# Ryan chose to promote it to the PRIMARY outcome definition (Option A). This cell re-runs, on that outcome, only the analyses that stay in the main text.
# It adds NO new kinds of analysis: every model below already exists in the notebook on the old definition. Analyses that do not use the baseline
# (absolute-rise criterion, dialysis, in-hospital death) are unchanged and are not repeated. Cells 13 and 14 are NOT needed.
# Aggregates only; counts below SMALL_CELL are masked. Paste back ONLY this cell's output.
import patsy
say("\n=== CELL 15: main-text analyses re-run on the HIERARCHY baseline outcome (post hoc) ===")
need = [c_ for c_ in ['aki_hier', 'base_hier', 'base_prior', 'peak_all', 'aki'] if c_ not in d.columns]
if need: raise SystemExit("CELL 15 needs Cell 12 to have run first in this session (missing: %s)." % need)
d['aki_h'] = d.aki_hier.astype(int); d['open'] = (d.approach == 'open').astype(int)
d['ratio_h'] = d.peak_all / d.base_hier
d['stage2plus_h'] = (d.ratio_h >= 2.0).astype(int); d['stage3_h'] = (d.ratio_h >= 3.0).astype(int)
say("CHECK: hierarchy AKI events", n_(d.aki_h.sum()), "(R11 and Cell 12: 227) | admissions with no baseline at all:", n_(d.base_hier.isna().sum()), "(R11: 0)")
say("Baseline source (rows = approach; columns 0 = first-24 h value, 1 = pre-admission median; small cells masked):"); say(cnt(pd.crosstab(d.approach, d.base_prior.notna().astype(int))).to_string())

def evalue15(r):
    r = 1/r if r < 1 else r
    return r + np.sqrt(r*(r-1))
def gcomp15(model, data):
    a = data.copy(); a['approach'] = 'open';         p1_ = model.predict(a).mean()
    b = data.copy(); b['approach'] = 'endovascular'; p0_ = model.predict(b).mean()
    return p1_, p0_

# ---- H1: the new primary ----
say("\n--- New primary: AKI by the hierarchy baseline ---")
m_h = fit(d, 'aki_h', T(), 'H1  PRIMARY (hierarchy baseline): modified Poisson RR, cluster by patient', 'poisson')
if m_h is None: raise SystemExit("The H1 model failed; send me the MODEL FAILED line above.")
est_h, lo_h, hi_h = RESULTS[-1]['est'], RESULTS[-1]['lo'], RESULTS[-1]['hi']
say("   (must equal P2 from Cell 12: 1.77 (1.41 to 2.22), 227 events)")
fit(d, 'aki_h', T(), 'H1b Same covariates, logistic odds ratio (secondary)', 'logit')
fit(d, 'aki_h', [T()[0]], 'H1u UNADJUSTED modified Poisson RR', 'poisson')
ols_ = smf.ols("aki_h ~ open", data=d).fit(cov_type='cluster', cov_kwds={'groups': d.subject_id.values})
say(f"H1v Unadjusted risk difference {100*ols_.params['open']:.1f} points ({100*ols_.conf_int().loc['open',0]:.1f} to {100*ols_.conf_int().loc['open',1]:.1f}), cluster-robust")
say("Crude AKI (hierarchy baseline): open", f"{d[d.approach=='open'].aki_h.mean():.1%}", "| endovascular", f"{d[d.approach=='endovascular'].aki_h.mean():.1%}", "| overall", f"{d.aki_h.mean():.1%}")
fml_h = "aki_h ~ " + " + ".join(T())
mb_h = smf.glm(fml_h, data=d, family=sm.families.Binomial()).fit()
q1, q0 = gcomp15(mb_h, d)
dd_ = d.reset_index(drop=True); grp_ = dd_.groupby('subject_id').indices; keys_ = np.array(list(grp_.keys())); rng_ = np.random.default_rng(20261027); bs_ = []
for _ in range(N_BOOT):
    idx_ = np.concatenate([grp_[k_] for k_ in rng_.choice(keys_, size=len(keys_), replace=True)]); s_ = dd_.iloc[idx_]
    try:
        a_, b_ = gcomp15(smf.glm(fml_h, data=s_, family=sm.families.Binomial()).fit(disp=0), s_); bs_.append((a_-b_, a_/b_))
    except Exception: pass
bs_ = np.array(bs_)
say(f"H1c Standardized (g-computation) risk if all open {q1:.1%} vs all endovascular {q0:.1%}; risk difference {100*(q1-q0):.1f} points "
    f"(95% CI {100*np.percentile(bs_[:,0],2.5):.1f} to {100*np.percentile(bs_[:,0],97.5):.1f}); standardized RR {q1/q0:.2f} ({np.percentile(bs_[:,1],2.5):.2f} to {np.percentile(bs_[:,1],97.5):.2f}); "
    f"{len(bs_)} patient-clustered bootstrap reps, seed 20261027")
lim_h = lo_h if est_h > 1 else hi_h
say(f"H1d E-value for RR {est_h:.2f}: {evalue15(est_h):.2f}; for the CI limit nearest the null ({lim_h:.2f}): {1.0 if (lo_h <= 1 <= hi_h) else evalue15(lim_h):.2f}")

# ---- severity on the hierarchy baseline ----
say("\n--- Severity (KDIGO ratio stages on the hierarchy baseline) ---")
sg_ = pd.cut(d.ratio_h, [0, 1.5, 2.0, 3.0, np.inf], right=False, labels=['no AKI', 'stage 1', 'stage 2', 'stage 3'])
say("Counts by approach (small cells masked):"); say(cnt(pd.crosstab(d.approach, sg_)).to_string())
fit(d, 'stage2plus_h', T(), 'H2  Stage 2 or 3 (ratio >= 2.0)', 'poisson')
fit(d, 'stage3_h', T(), 'H3  Stage 3 (ratio >= 3.0)', 'poisson')

# ---- the sensitivity analyses that stay in the main text ----
say("\n--- Main-text sensitivity analyses, repeated on the hierarchy outcome ---")
first_ = d.sort_values('admittime').drop_duplicates('subject_id')
fit(first_, 'aki_h', T(), 'H4  First admission per patient only', 'poisson')
fit(d[d.rupture == 0], 'aki_h', [t_ for t_ in T() if t_ != 'rupture'], 'H5  Excluding ruptured repairs', 'poisson')
fit(d, 'aki_h', T(['vasopressor']), 'H6  Add any-vasopressor (post-procedure variable, mediator risk)', 'poisson')
if 'log_cr_per_day' in d.columns: fit(d, 'aki_h', T(['log_cr_per_day']), 'H7  Add log(creatinine draws per day)', 'poisson')
else: say("H7 skipped: log_cr_per_day not found (it is created in Cell 2).")
def ow_y(df, ycol):
    Xd = patsy.dmatrix(COV_FORMULA, df, return_type='dataframe')
    ps_ = sm.GLM(df['open'], Xd, family=sm.families.Binomial()).fit(disp=0).predict(Xd).values
    w_ = np.where(df['open'] == 1, 1-ps_, ps_); y_ = df[ycol].values; o_ = df['open'].values
    return (w_*y_*o_).sum()/(w_*o_).sum(), (w_*y_*(1-o_)).sum()/(w_*(1-o_)).sum()
r1_, r0_ = ow_y(d, 'aki_h'); rng2_ = np.random.default_rng(20261028); bo_ = []
for _ in range(500):
    idx_ = np.concatenate([grp_[k_] for k_ in rng2_.choice(keys_, size=len(keys_), replace=True)])
    try: a_, b_ = ow_y(dd_.iloc[idx_], 'aki_h'); bo_.append((a_/b_, a_-b_))
    except Exception: pass
bo_ = np.array(bo_)
say(f"H8  Overlap-weighted risks: open {r1_:.1%}, endovascular {r0_:.1%}; RR {r1_/r0_:.2f} ({np.percentile(bo_[:,0],2.5):.2f} to {np.percentile(bo_[:,0],97.5):.2f}); "
    f"RD {100*(r1_-r0_):.1f} points ({100*np.percentile(bo_[:,1],2.5):.1f} to {100*np.percentile(bo_[:,1],97.5):.1f}); 500 patient-clustered bootstrap reps")
for L_ in ['abdominal', 'thoracic', 'thoracoabdominal']:
    sub_ = d[d.location == L_]
    fit(sub_, 'aki_h', [t_ for t_ in T() if 'location' not in t_], f'H9  Stratum: {L_} (N={n_(len(sub_))})', 'poisson')
fi_ = "aki_h ~ C(approach, Treatment('endovascular'))*C(location, Treatment('abdominal')) + age_at_admission + C(gender) + cci + C(anchor_year_group) + rupture + C(urg, Treatment('elective'))"
try:
    m0_ = smf.glm(fi_, data=d, family=sm.families.Poisson()).fit()
    mi_ = smf.glm(fi_, data=d, family=sm.families.Poisson()).fit(cov_type='cluster', cov_kwds={'groups': d.loc[m0_.model.data.row_labels, 'subject_id'].values})
    its_ = [k_ for k_ in mi_.params.index if ':' in k_]
    R_ = np.zeros((len(its_), len(mi_.params)))
    for i_, k_ in enumerate(its_): R_[i_, list(mi_.params.index).index(k_)] = 1
    say(f"H10 Joint Wald test of approach x location: p={float(mi_.wald_test(R_, scalar=True).pvalue):.3f}")
except Exception as e_:
    say("H10 interaction model failed:", type(e_).__name__, str(e_)[:80])

say("\nREFERENCE (old definition, ledger): primary RR 1.87 (1.49 to 2.34), 229 events; logistic OR 2.94 (1.96 to 4.43); unadjusted RR 1.70 (1.37 to 2.10); unadjusted RD 17.3;")
say("standardized RD 21.13 (13.39 to 28.62); E-value 3.14 and 2.34; stage 2 or 3 RR 1.93 (1.37 to 2.71), 131 events; stage 3 RR 0.97 (0.55 to 1.71), 66 events;")
say("first admission 1.91 (1.51 to 2.41); excluding rupture 1.80 (1.36 to 2.37); vasopressor 1.62 (1.29 to 2.03); draws per day 1.71 (1.38 to 2.12); overlap 1.85 (1.51 to 2.32);")
say("strata abdominal 1.70 (1.29 to 2.23), thoracic 2.77 (1.60 to 4.81), thoracoabdominal 1.80 (0.97 to 3.34); interaction joint p=0.905.")
out15 = pd.DataFrame([r for r in RESULTS if str(r['model']).startswith('H')])
if len(out15):
    out15['Est (95% CI)'] = out15.apply(lambda r: f"{r.est:.2f} ({r.lo:.2f} to {r.hi:.2f})", axis=1)
    out15['p'] = out15.p.apply(lambda x: '' if pd.isna(x) else ('<0.001' if x < 0.001 else f"{x:.3f}")); out15['events'] = out15.events.apply(lambda v: n_(v))
    print("\n=================== CELL 15 RESULTS TABLE ==================="); print(out15[['model', 'N', 'events', 'measure', 'Est (95% CI)', 'p']].to_string(index=False)); print("=================== END ===================")
