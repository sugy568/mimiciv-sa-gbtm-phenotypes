# mimiciv-sa-gbtm-phenotypes

Analysis code for the manuscript:

> **Joint lactate–creatinine trajectories identify three phenotypes with markedly different mortality in *Staphylococcus aureus* infection: a retrospective cohort study using MIMIC-IV**

A retrospective cohort study based on MIMIC-IV v3.1 (2008–2022). Adult admissions with a
first microbiological detection of *S. aureus* (n = 13,773) were phenotyped by group-based
trajectory modelling (GBTM, `lcmm::hlme`) of serial lactate and creatinine values measured
within 0–72 h of the first detection. A joint three-level phenotype (P1 low → P3 high) was
derived by summing the severity ranks of the two per-marker class assignments. Outcomes were
28-day mortality (Cox models), competing-risk cumulative incidence (Aalen–Johansen), and the
incremental discrimination of the phenotype over clinical baselines and single-point
biomarkers (ΔC with DeLong tests, IDI, NRI).

## Data availability and DUA compliance

The source data are from **MIMIC-IV version 3.1**, publicly available to credentialed users
through PhysioNet: <https://physionet.org/content/mimiciv/3.1/> (dataset DOI
[10.13026/kpb9-mt58](https://doi.org/10.13026/kpb9-mt58)). Access requires registration,
signature of a data use agreement (DUA), and completion of training in the responsible
conduct of research (CITI Program).

> **No patient-level data are included in this repository.** The PhysioNet DUA prohibits
> redistribution of patient-level data. Users must obtain their own credentialed access and
> regenerate the analysis tables with the SQL scripts in [`sql/`](sql/). Only aggregate
> summary tables, trace files and figures are committed here.

## Environment

| Component | Version used |
|---|---|
| R | 4.3.3 (2024-02-29) |
| lcmm (GBTM fitting) | 2.2.2 |
| pROC (ROC / DeLong) | 1.19.1 |
| survival (Cox models) | 3.5-8 |
| cmprsk (competing risks) | 2.2-12 |
| ggplot2 / gridExtra | 4.0.3 / 2.3.1 |
| PostgreSQL + MIMIC-IV | 12+ with `mimiciv_hosp` / `mimiciv_icu` / `mimiciv_derived` loaded |

Full dependency lockfile: [`package_versions.csv`](package_versions.csv) and
[`sessionInfo.txt`](sessionInfo.txt). Missing packages can be installed with
`Rscript R/00_install_packages.R`.

## Repository structure

```
mimiciv-sa-gbtm-phenotypes/
├── sql/            # MIMIC-IV extraction scripts (PostgreSQL), run in numeric order
│   ├── 01_sa_cohort_feasibility.sql
│   ├── 02_sa_specimen_resistance_medications.sql
│   ├── 03_gbtm_trajectory_dataset.sql
│   ├── 04_gbtm_marker_intersection_fast.sql
│   ├── 05_sa_antibiotics_gsn.sql
│   ├── 06_optional_emar_administration.sql
│   └── deprecated/                     # earlier-generation scripts, kept for provenance
├── R/              # analysis pipeline (run from this folder)
│   ├── 00_install_packages.R           # dependency installation
│   ├── 01_load_data.R                  # semantic, filename-agnostic data loader
│   ├── 00_RUN_ALL.R                    # one-key master: GBTM refit → all tables & figures
│   ├── 20_manuscript_tables_figures.R  # all manuscript tables and figures
│   ├── 21_rrt_sensitivity.R            # RRT-definition sensitivity analysis
│   └── 40_incremental_C_statistic.R    # incremental discrimination (ΔC, IDI, NRI)
├── data/           # ← place the four exported CSVs here (git-ignored, never commit)
├── figures/        # manuscript figures Fig01–08 and FigS01–S11 (PNG/PDF/SVG)
├── tables/         # manuscript table outputs (aggregate only) + incremental_C/
└── trace/          # model-selection trace, input fingerprints, key numbers
```

## Reproduction

### Step 1 — SQL extraction (PostgreSQL)

Execute the scripts in [`sql/`](sql/) **in numeric order** against a PostgreSQL instance
with MIMIC-IV v3.1 loaded. Scripts are idempotent and print expected intermediate row
counts so partial runs can be verified. Export the four analysis tables (cohort master,
0–72 h marker long table, susceptibility long table, RRT indicators) to `data/` as UTF-8
CSVs. Do **not** commit these files and avoid non-ASCII paths.

### Step 2 — Main pipeline

```bash
Rscript R/01_load_data.R               # sanity check of file location and columns
Rscript R/00_RUN_ALL.R                 # one-key run, ~12–20 min (excl. GBTM grid search)
Rscript R/00_RUN_ALL.R --stage=gbtm    # refit phenotypes only
Rscript R/00_RUN_ALL.R --stage=plot    # tables/figures only (reuse phenotype file)
```

GBTM specification: `hlme(value ~ time_s + I(time_s^2), random = ~time_s,
mixture = ~time_s + I(time_s^2), subject = "hadm_id", ng = K)` with `time_s = time_h / 72`;
K = 1–4 screened by BIC with a 30-start grid search; every fitting attempt is logged to
`trace/gbtm_search_log.csv`.

### Step 3 — Incremental discrimination analysis

```bash
Rscript R/40_incremental_C_statistic.R
```

Self-contained script; all settings are environment-variable overridable (`DATA_DIR`,
`OUT_DIR`, `PHENO_MODE`, `GRID_REP`, `SEED`, `BOOT_B`). Outputs nested-model C-statistic
comparisons (ΔC + bootstrap CI + DeLong test + LR χ² + IDI/NRI + Brier + AIC) for four
pre-specified comparator scenarios.

Outputs of both steps are written to `out_manuscript/{figures,tables,trace}`.

## Notes for reviewers

- `trace/gbtm_search_log.csv` — every GBTM fitting attempt (K, seeds, BIC, convergence),
  documenting the model-selection process.
- `trace/run_manifest.txt` — MD5 fingerprint and dimensions of the exact input file used
  for the manuscript run; `trace/key_numbers.txt` — key result numbers.
- Committed outputs under `figures/`, `tables/` and `trace/` correspond to the manuscript run.

## License

Code is released under the MIT License (see [LICENSE](LICENSE)).
MIMIC-IV data are **not** included and are not redistributable; access them via PhysioNet
under your own data use agreement.
