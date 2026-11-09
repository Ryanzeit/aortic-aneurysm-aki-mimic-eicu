# Open versus endovascular aortic aneurysm repair and acute kidney injury in intensive care

Analysis code for the manuscript *Open versus endovascular aortic aneurysm repair and acute kidney injury in intensive care: a retrospective cohort study using MIMIC-IV, with an exploratory eICU-CRD analysis*.

The repository holds **code only**. It contains no patient-level data, no data extracts and no notebook outputs.

## Data access

The data are not included and cannot be shared here. Both databases are credentialed PhysioNet resources:

- MIMIC-IV (hospital and ICU modules): https://physionet.org/content/mimiciv/
- eICU Collaborative Research Database (eICU-CRD): https://physionet.org/content/eicu-crd/

To run the code you must hold your own PhysioNet credentialing, complete the required training, and sign the data use agreements. Do not upload the data, or any file derived from it at patient level, to this repository (`.gitignore` blocks common data file types).

## How the code is organised

The work was developed interactively in Google Colab. Paths follow the Colab layout (`/content/drive/MyDrive/mimic-iv/...`, `/content/drive/MyDrive/eICU-CRD/...`); edit the path variables at the top of each file to match your own copy of the data. All notebooks are saved with outputs cleared. Printed results in the scripts are aggregates, with counts below 11 masked as `<11` (a conservative precaution set by the authors, not a PhysioNet rule).

| Folder | Contents |
|---|---|
| `01_cohort_build/` | Cohort construction for both databases, in run order: `01_mimic_cohort_build` (MIMIC-IV procedure, ICU stay, creatinine, comorbidity and vasopressor data), `02_mimic_exclusions_eicu_merge_and_checks` (rupture and urgency, ESRD exclusion, no-creatinine exclusion, pooled interaction check, Table 1, supplement code lists), `03_eicu_cohort_build` (eICU-CRD cohort, location, age, AKI). |
| `02_mimic_analysis/` | Primary and sensitivity analyses in MIMIC-IV. `AKI_notebook5_v2_6_cells_0-11` first (use Run all), then the add-on cells in the same session: `AKI_notebook5_v2_7_addon_cells_12-14` (Cell 12 builds the pre-admission baseline hierarchy that is the primary definition), `cell15_rerun_on_hierarchy_baseline.py` (main-text models on the primary definition) and `cell16_creatinine_range_check.py` (descriptive range and duplicate check). |
| `03_eicu_analysis/` | `eicu_step2_location_model.py`: thoracic versus abdominal location and AKI in eICU-CRD with the text-derived comorbidity count. |
| `04_verification/` | Independent recomputation of the main results with a different code path (DuckDB SQL plus a hand-written Poisson with patient-clustered standard errors): `R10` (first definition), `R11` (pre-admission baseline), `R12` (primary definition), and `cell12_verification_checks.py` (handling of the one admission without a baseline, recount of the pre-procedure baseline subset). |
| `05_figures/` | Scripts for Figure 1 (flow diagram) and Figure 2 (forest plot). |

### Run order

1. `01_cohort_build/01` then `02`; `03` for eICU-CRD (the notebooks write intermediate cohort and creatinine files that later steps read; file names are in the code).
2. `02_mimic_analysis/AKI_notebook5_v2_6_cells_0-11.ipynb` (Run all), then Cells 12 to 14, `cell15_...` and `cell16_...` in the same session. These cells reuse variables from the earlier cells (`d`, `fit`, `say`, `T`).
3. `03_eicu_analysis/eicu_step2_location_model.py`.
4. `04_verification/` scripts (each is self-contained).
5. `05_figures/` scripts (hard-code the published estimates; no data needed).

### Notes

- Cell numbering in comments follows the interactive analysis session, so comments that mention a cell or a version refer to that session. The verification file `cell12_verification_checks.py` is separate from the add-on notebook's Cell 12 (baseline lookback).
- The primary outcome is KDIGO creatinine-ratio AKI (peak at least 1.5 times baseline), with baseline = median creatinine 7 to 365 days before admission where available, otherwise the lowest value in the first 24 hours. The first-24-hour baseline is a sensitivity analysis.
- Analyses are associational and post hoc; there was no registered protocol.

## Requirements

Python 3 with the packages in `requirements.txt` (Colab defaults plus `duckdb` and `comorbidipy`).

## Citation

See the manuscript for the full citation. A Zenodo archive of a tagged release provides a permanent DOI for this code.

## Licence

MIT (see `LICENSE`).
