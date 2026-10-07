# Kora calibration (2025-11-01)

> **Note:** This is research code that is functional but does not follow best practices.
> The structure is documented here to help readers reproduce the analysis as well as
> possible. My apologies to the reader who has to figure it out, as I do not have time to
> clean it all up. If it is any consolation, I hate the state of the code very much.

Future efforts will incorporate the ensemble calibration and assessment process into
[Kora.jl](https://github.com/open-AIMS/Kora.jl) or a companion package.

---

## Purpose

This study calibrates and evaluates the [Kora.jl](https://github.com/open-AIMS/Kora.jl)
individual-based coral reef model against field observations for two reefs:

- **Moore Reef** (16071S): Offshore North region, Great Barrier Reef
- **Masig Reef**: Torres Strait region

The workflow first processes raw EcoRRAP benthic survey data into annual mean cover estimates by
functional group, then performs an exploratory sensitivity analysis of the growth and survival
functions, fits reef-specific regression models to the field data, calibrates an ensemble of
plausible initial reef states against historical observations, and finally performs sensitivity
analysis to identify which parameters most influence model behaviour.

### Unified dataset (all sites)

The all-sites collation workflow has been moved to a separate project
(`EcoRRAP-data-collator`). It produces `ecorrap_unified.parquet`, which is consumed by the
Stage 1–2 scripts in this repository.

---

## Requirements

### Julia

Julia ≥ 1.12 is recommended. Install dependencies by activating the project and running:

```julia
] instantiate
```

This will install all required dependencies.

### Key packages (declared in `common.jl`)

| Package | Purpose |
|---|---|
| `Kora` | Core reef simulation model (separate repository) |
| `BlackBoxOptim` | Evolutionary optimisation for calibration |
| `QuasiMonteCarlo` | Sobol sampling for sensitivity analysis |
| `HypothesisTests` | KS-test used internally by the PAWN implementation |
| `PairPlots` | Parameter correlation visualisation |
| `CairoMakie` | Figures |
| `Distributed` | Parallel calibration workers |

### Parallel execution

Scripts `03a` and `04a` spawn worker processes via `Distributed.addprocs`. The number of
workers is set at the top of each script (`n_workers = 20`). Adjust this to match available CPU
cores before running.

---

## External data

The following files are **not** included in this repository and must be obtained separately before
running any scripts.

### Restricted data (not publicly available)

These files were provided under data sharing arrangements and cannot be redistributed. They
are expected in the study `data/` directory.

| File | Source |
|---|---|
| `ecorrap_adult_juv_combined_2021_2023_24062025.csv` | EcoRRAP program (AIMS) |
| `EcoRRAP data for IPM_250624.csv` | EcoRRAP program (AIMS) |
| `ecorrap_to_groups.csv` | EcoRRAP program (AIMS) - maps taxon codes to Kora functional groups |

The first incorporates juvenile quadrat data from:
> Doropoulos, C., Alvarez-Noriega, M., Fabricius, K., Ferrari, R., Mumby, P.J., Noonan, S.H.C., Orr, M., >
> Salee, K., 2025. Impact of environmental gradients on juvenile coral demography across the Great Barrier
> Reef and Torres Strait. Coral Reefs. https://doi.org/10.1007/s00338-025-02742-6

The second dataset was collated by Dr. Anna Cresswell (Australian Institute of Marine Science)
and can be provided on request.

The third file was compiled based on advise from Dr. Anna Cresswell and Dr. Renata Ferrari.

#### Processed outputs derived from restricted data

The following files - produced by stages 1 and 2 from the restricted inputs above - are
included in this repository. **Stages 3 - 5 (calibration and assessment) can therefore be
run without access to the raw restricted data.**

| Location | Contents |
|---|---|
| `data/offshore_north/overall/` | `offshore_north_growth.csv`, `offshore_north_survival.csv`, `*_growth_models.json`, `*_survival_models.json` |
| `data/offshore_north/moore/` | `offshore_north_moore_growth.csv`, `offshore_north_moore_survival.csv`, `*_growth_models.json`, `*_survival_models.json` |
| `data/torres_strait/overall/` | `torres_strait_growth.csv`, `torres_strait_survival.csv`, `*_growth_models.json`, `*_survival_models.json` |
| `data/torres_strait/masig/` | `torres_strait_masig_growth.csv`, `torres_strait_masig_survival.csv`, `*_growth_models.json`, `*_survival_models.json` |

The `*_models.json` files are the fitted growth and survival regression models consumed by
the calibration scripts. The `*_growth.csv` / `*_survival.csv` files are the pre-processed
individual-level training data used to fit them.

The [Canonical Reefs](https://github.com/gbrrestoration/canonical-reefs) dataset is used to
identify reefs of interest by their UNIQUE IDs.

Manta Tow data were sourced directly from the
[AIMS Reef Monitoring Dashboard](https://apps.aims.gov.au/reef-monitoring/reefs) for the
reefs of interest.

#### Unified dataset — additional restricted files

The following additional files are required for the all-sites collation workflow in
`EcoRRAP-data-collator`. They are also placed in the study `data/` directory.

| Path (relative to study root) | Contents |
|---|---|
| `data/ecorrap_benthic/latest/cover_estimate_DESCRIPTION.csv` | Photoquadrat benthic cover (%) per species and transect for all 13 sites |
| `data/ecorrap_benthic/latest/EcoRRAP_Community_Composition_Labelset.xlsx` | Updated labelset mapping DESCRIPTION codes to functional groups; export to CSV and update the collator config if the current CSV is missing species |
| `data/ecorrap_benthic/EcoRRAP_Labelset_mapping_Maren_Toor_2026-01-20.csv` | Older labelset mapping used by the collator as the default; one-to-one DESCRIPTION → functional group, no ambiguous entries |

The collator also reads EcoRRAP oceanographic logger data (temperature, salinity, current
speed, depth, waves and PAR). These files are kept with `EcoRRAP-data-collator` and not in this
repository. Source: Australian Institute of Marine Science (AIMS) (2023), *EcoRRAP oceanographic
logger data*, https://apps.aims.gov.au/metadata/view/778cc4e3-7979-4601-8fdf-772f345d8956
(CC BY-NC 3.0 AU).

The site code lookup table (`reef_site_code_lookup.csv`) is generated by
`EcoRRAP-data-collator` and is not included in this repository. It maps reef names to the site
codes used in the benthic and oceanographic datasets (which use different codes for some
Southern GBR sites).

### Publicly available data

These files must be placed in the locations shown relative to the study folder.

| File | Location | Source |
|---|---|---|
| `Moore Reef_Manta Tow_line_chart_modelled_2025-12-28.csv` | `data/` | AIMS LTMP |
| `dhw_scens.nc` | `data/dhw/` | RRAP / CoralBlox scenario data |
| `DHW_Moore_Reef.csv`, `DHW_Masig_reef.csv` | `data/dhw/` | NOAA CRW or extracted from above |
| `rrap_canonical_2025-07-15-T10-48-29.gpkg` | `data/` | RRAP canonical reef dataset |

`data/dhw/DHW_KBN_extracted.csv` holds the combined Moore and Masig DHW series from which the
two per-reef files were split. No script reads it; it is kept for provenance.

---

## Workflow

Scripts are in the `scripts/` subdirectory and numbered in execution order.
**Run all scripts from within `scripts/`:**

If running directly from the command line:

```
cd scripts
julia --project=.. 00a_prep_ecorrap_data.jl
```

or via the REPL (from the project root):

```
; cd scripts
include("00a_prep_ecorrap_data.jl")
```

Within each stage, `a` scripts produce data consumed by `b` scripts.

To run the pipeline through Stage 5 in order, use `julia --project=.. run_all.jl`
from `scripts/`. The runner logs failures and continues to later stages, so check its
summary for failed scripts.

### Stage 0: Data preparation

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/00a_prep_ecorrap_data.jl` | `data/ecorrap_benthic/latest/cover_estimate_DESCRIPTION.csv`, labelset mapping | `data/ecorrap_benthic/moore_estimate.csv`, `masig_estimate.csv`, `data/Masig_Reef_EcoRRAP_estimate.csv` | Produces annual mean cover by functional group for Moore and Masig from the expanded benthic data. Original version (old data format) archived as `00a_prep_ecorrap_data_ARCHIVED.jl`. Required for Stages 1–5. |
| `scripts/00b_study_area_map.jl` | `data/rrap_canonical_2025-07-15-T10-48-29.gpkg`, `data/EcoRRAP data for IPM_250624.csv`, `data/ecorrap_unified.parquet` | `figs/study_area_overview.png` | Produces the GBR + Torres Strait study context figure with regional insets and highlighted study reefs. |

### Stage 1: Exploratory sensitivity analysis of growth/survival functions

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/01a_exploratory_SA.jl` | `data/ecorrap_unified.parquet` + species map CSV | `figs/sensitivity/` PAWN heatmaps | Region-wide SA; helper functions used by `01b`/`01c` |
| `scripts/01b_exploratory_SA_offshore_north.jl` | `data/ecorrap_unified.parquet` + species map CSV | `data/offshore_north/{overall,moore}/` CSVs + model `.json` files, `figs/sensitivity/offshore_north/` | Fits and evaluates growth/survival models for offshore north region and Moore Reef |
| `scripts/01c_exploratory_SA_torres_strait.jl` | `data/ecorrap_unified.parquet` + species map CSV | `data/torres_strait/{overall,masig,...}/` CSVs + model `.json` files, `figs/sensitivity/torres_strait/` | Same for Torres Strait |
| `scripts/01d_compile_SA_figures.jl` | Sensitivity PNGs from `01a`-`01c` | `figs/sensitivity/combined/SA_*.png` | Combines growth and survival plots into four paper figures |
| `scripts/01e_paper_dataset_summary.jl` | Region-level growth and survival CSVs from `01b`/`01c` | `data/logger_report/paper_dataset_summary.md`, `paper_table1_dataset_overview.md` | Generates dataset diagnostics and the Table 1 body included by the paper |

### Stage 2: Regression model fitting (diameter-based)

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/02a_offshore_north_fit_to_diameter.jl` | `data/ecorrap_unified.parquet` + species map CSV | `data/offshore_north/{overall,moore}/*_models.json`, `figs/regressions/offshore_north/` | Fits final growth (degree-1) and survival (degree-2) polynomial regressions; **these `.json` files are the models used in calibration** |
| `scripts/02b_torres_strait_fit_to_diameter.jl` | `data/ecorrap_unified.parquet` + species map CSV | `data/torres_strait/{overall,masig}/*_models.json`, `figs/regressions/torres_strait/` | Same for Torres Strait |

### Stage 3: Moore Reef ensemble calibration and assessment

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/03a_moore_ensemble.jl` | Config `scripts/_16071S_config.jl`, model `.json` files, LTMP manta tow CSV, DHW CSV | `data/ensemble/offshore_north/moore/16071S_*` (`.dat`, `.h5`, `.csv`), calibration figures in `figs/`, summary CSVs | Runs parallel multi-start ensemble calibration; see [Iterative calibration](#iterative-calibration) |
| `scripts/03b_moore_ensemble_assessment.jl` | Output of `03a`, LTMP data, EcoRRAP benthic estimate | `data/sensitivity/offshore_north/moore/ensemble/16071S_*` (`.h5`, `.csv`), `data/ensemble/offshore_north/moore/16071S_parameter_correlations.csv`, sensitivity figures | Unconstrained and ensemble-constrained PAWN SA; temporal and lagged sensitivity analysis |

### Stage 4: Masig Reef ensemble calibration and assessment

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/04a_masig_torres_strait_ensemble.jl` | Config `scripts/_masig_config.jl`, model `.json` files, EcoRRAP benthic CSV, DHW CSV | `data/ensemble/torres_strait/masig/masig_*` (`.dat`, `.h5`, `.csv`), calibration figures, summary CSVs | Same parallel calibration workflow as `03a` |
| `scripts/04b_masig_ensemble_assessment.jl` | Output of `04a`, EcoRRAP benthic CSV | `data/sensitivity/torres_strait/masig/ensemble/masig_*` (`.h5`, `.csv`), `data/ensemble/torres_strait/masig/masig_parameter_correlations.csv`, sensitivity figures | Same SA workflow as `03b` |

### Stage 5: Combined sensitivity heatmaps (Moore + Masig)

| Script | Inputs | Key outputs | Notes |
|---|---|---|---|
| `scripts/05_combined_sa_heatmaps.jl` | Cached PAWN results from `03b` (`16071S_*_pawn_results.h5`) and `04b` (`masig_*_pawn_results.h5`) | `figs/sensitivity/combined/` | Combines Moore and Masig PAWN heatmaps into shared figures on a common colour scale; no calibration, plotting only |

---

## Configuration

Reef-specific settings for Moore and Masig are kept in standalone config files that are `include`d
at the top of each calibration script:

| File | Reef | Controls |
|---|---|---|
| `scripts/_16071S_config.jl` | Moore Reef | `ReefConfig` (area, depth, density, initial proportions, excluded years, disturbance years), `CalibrationDataPaths`, `OptimizationConfig`, `CalibrationSettings`, `param_bounds` |
| `scripts/_masig_config.jl` | Masig Reef | Same structure |

**To adapt to a new reef**, copy one of the config files, update all fields, and create new `a`/`b`
script pair following the existing pattern.

### Key configuration parameters

| Parameter | Location | Effect |
|---|---|---|
| `exclude_years` | `ReefConfig` | Years removed from calibration (held out for validation or disturbed) |
| `disturbance_years` | `ReefConfig` | Marked on time-series plots as known disturbance events |
| `fitness_threshold` | `OptimizationConfig` | Scores below this are accepted as ensemble members; lower = stricter |
| `ensemble_members` | `OptimizationConfig` | Target number of accepted parameter sets |
| `max_steps` | `OptimizationConfig` | Maximum optimisation steps per trial |
| `scalers_*` bounds | `param_bounds` | Range allowed for per-group growth rate multipliers; narrow these (e.g. `(0.95, 1.05)`) if post-calibration growth is unrealistic |

---

## Iterative (multi-start ensemble) calibration

The `03a`/`04a` scripts implement a deliberate multi-round workflow:

1. **First run** — no prior candidates exist; the script runs from scratch and saves
   `*_initial_guess.h5` alongside the ensemble output.
2. **Subsequent runs** — the script detects `*_initial_guess.h5` and seeds the initial
   population with previously found candidates, improving convergence.
3. **Resuming** — if `*_optim_state.dat`, `*_optim_best.h5`, and `*_tracked_candidates.h5`
   all exist, the script loads them directly and skips re-optimisation, proceeding straight to
   visualisation.

To force a fresh run, delete or rename the cached `.dat` and `.h5` files in the relevant
`data/ensemble/<region>/<reef>/` directory.

Future improvements will aim to automate the process.

---

## Output directory layout

```
scripts/    # all runnable .jl scripts and config files
data/
├── <region>/<reef_or_overall>/     # regression model .json files and fitted-data .csv files
├── ecorrap_benthic/                # EcoRRAP benthic cover data and estimates
│   ├── moore_estimate.csv          # annual mean cover by functional group — Moore Reef
│   ├── masig_estimate.csv          # annual mean cover by functional group — Masig Reef
│   └── latest/                     # expanded benthic cover input files (all 13 sites)
├── ecorrap_unified.parquet         # unified IPM + benthic cover + oceanographic dataset (produced externally)
├── dhw/                            # degree heating week time series
├── ensemble/<region>/<reef>/       # calibration outputs (*_ensemble_output.dat, *_initial_guess.h5, etc.)
└── sensitivity/<region>/<reef>/ensemble/   # sensitivity analysis samples and PAWN results

figs/
├── regressions/<region>/<reef_or_overall>/   # growth and survival regression diagnostic plots
├── sensitivity/<region>/<reef_or_overall>/   # exploratory PAWN sensitivity figures
│   └── ensemble/                             # ensemble-constrained SA figures
└── <reef_id>_calib_param_interactions.png    # parameter–fitness scatter plots
```

---

## Non-obvious design choices

- **PAWN sensitivity analysis is a custom implementation** in `src/sensitivity.jl`, ported from
  [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl), which itself was adapted from the
  [SALib Python package](https://salib.readthedocs.io). It is not sourced from any Julia
  sensitivity analysis package.

- **Group proportions are not parameterised directly.** The five proportion parameters are Gamma
  quantiles that are transformed via an internal `gamma_to_dirichlet` function to ensure they sum
  to 1. Do not interpret their raw values as proportions.

- **The fitness metric is a composite score** equal to NKGE (normalised bias β, normalised
  variability ratio α, and Pearson correlation term) plus an endpoint MAE penalty, a low-cover
  penalty, and a benthic rank score. Lower is better. The exact formulation is in
  `src/calibration_helpers.jl`.

- **`n_workers = 20`** at the top of parallel scripts is hardware-specific. Set it to the number
  of physical cores available, minus one or two for the OS.

- **Growth scalers for Masig are fixed at effectively 1.0** (bounds `(0.999, 1.001)` in
  `scripts/_masig_config.jl`) and so are not calibrated. The EcoRRAP-derived growth functions
  already encode reef-specific rates, and letting the scalers move allows the optimiser to suppress
  growth during the calibration period and then produce unrealistic acceleration once DHW stress
  drops. The bounds cannot be an exact point because uniform sampling and the nearest-neighbour
  diversity normalisation both need a nonzero width. Moore's scalers are calibrated over wider,
  per-group ranges (see `scripts/_16071S_config.jl`).
