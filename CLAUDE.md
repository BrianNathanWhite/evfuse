# evfuse

Spatial data fusion for extreme value analysis. Implements a two-stage frequentist framework for fusing sparse observations and dense simulations, developed for U.S. coastal sea levels (White, Blanton, Luettich, & Smith, in preparation) but applicable to any annual maxima from multiple spatial data sources. Extends the Russell et al. (2019) coregionalization approach with partial observations and analytic gradients.

## Key conventions
- Parameter-major ordering throughout (parameter index varies fastest)
- Log-rho parameterization for optimization (ensures rho > 0)
- Lambda = 300 km Wendland C2 taper for bootstrap W
- B = 500 bootstrap replicates
- 20 multi-start random restarts from log-uniform[50, 5000] km
- Prediction grid: filtered coastal points at 15 km spacing, lon -98 to -66, lat 24 to 46
- Seeds: bootstrap (42), retry loop (101-105), multi-start (2026), return level simulation (123)
- `source_params` abstraction: named list mapping source labels to parameter indices (e.g., `list(NOAA = 1:3, ADCIRC = 4:6)`); all core functions accept this via `dat$source_params`

## Project structure

### R source files
- R/data.R — data loading (`load_data()` with `source_params`), distance matrix, `gulf_data` docs
- R/gev.R — Stage 1 pointwise GEV fitting via extRemes::fevd
- R/covariance.R — exponential covariance, Wendland taper, Kronecker construction
- R/bootstrap.R — block bootstrap for W, tapering (`taper_W`), embedding (`embed_W`)
- R/spatial_model.R — Stage 2 MLE with analytic gradient (`fit_spatial_model`)
- R/naive_model.R — 3-dim baseline models (NOAA-only, ADCIRC-only)
- R/kriging.R — GP kriging with partial observations (`predict_krig`)
- R/return_levels.R — 100-year return levels via delta method and simulation
- R/loo_cv.R — closed-form block LOO-CV (Rasmussen & Williams 2006)
- R/trend_diagnostics.R — Mann-Kendall test, GEV trend test
- R/detrend.R — nonstationary GEV for detrending robustness check

### Scripts
- scripts/run_pipeline.R — full pipeline (Stage 1 -> bootstrap -> Stage 2 -> krige -> return levels)
- scripts/run_naive.R — fit and evaluate naive baseline models
- scripts/run_loo.R — LOO-CV comparison (joint vs NOAA-only)
- scripts/run_detrended.R — detrended pipeline robustness check
- scripts/run_trends.R — trend diagnostics at NOAA and ADCIRC sites
- scripts/multi_start.R — multi-start optimization for Stage 2
- scripts/plot_results.R — main visualization (Stage 1 maps, kriged params, return levels)
- scripts/plot_comparison.R — comparison plots and tables (joint vs NOAA-only vs ADCIRC-only)
- scripts/print_tables.R — formatted console output of comparison tables
- scripts/save_tables.R — save comparison tables as RDS
- scripts/regression_check.R — QC regression check (18 tests against published numbers)

### Package infrastructure
- data/gulf_data.rda — bundled 129-site example dataset
- README.md — package overview, installation, quick start, customization
- vignettes/evfuse-tutorial.Rmd — tutorial vignette (pre-computed, eval=FALSE)
- tables/table_detrend_robustness.tex — LaTeX supplementary table

### Other
- data-raw/ — source data (.csv), fitted models (.rds), prediction grid, diagnostics
- figures/ — all output plots (13 PNGs)
- tests/testthat/test-core.R — 33 unit tests
- evfuse_notes.tex — running working notes with analytical descriptions of all results

## Data
- Source: `data-raw/combined_max.csv` or `data(gulf_data)`
- 129 sites: 29 NOAA, 100 ADCIRC
- 34-43 years of annual maxima per site

## Development practices
- Use `devtools::load_all()` and `devtools::test()` to iterate
- Run `test()` after every change
- Don't install packages without asking first
- If a fix involves changing mathematical structure (not just R syntax), explain what's changing and why before implementing
- Don't add co-authored-by lines to commits

## Implementation notes
- `extRemes::fevd` returns observed information matrix directly (positive definite), not Hessian of negative log-likelihood. Use `solve(fit$results$hessian)` directly for vcov.
- Delta method for log_sigma transformation: J = diag(1, 1/sigma, 1)
- Rho optimized in log space: `pack_params` stores `log(rho)`, `unpack_params` exponentiates. Gradient uses chain rule.
- Analytic gradient contracts Q = V_inv - alpha*alpha^T with Omega_i into pxp H matrices. See `data-raw/evfuse_gradient.pdf`.
- Sigma_obs built directly in observed space (never form full Lp x Lp).
- nll/nll_grad share a cache keyed on parameter vector to avoid redundant Cholesky.
- All core functions generalized via `source_params`; defaults match NOAA/ADCIRC case.

## Current state (Feb 2026)

### Completed
- Stage 1: All 129/129 sites converge, all vcov matrices PSD
- Stage 2: Joint 6-dim model fitted (NLL = -103.61, 33 params)
- Naive baselines: NOAA-only (NLL = +13.67) and ADCIRC-only (NLL = -75.49)
- Kriging + 100-year return levels for all three models
- LOO-CV: Joint wins 28/29 NOAA sites, LPD gain +32.4, RL RMSE reduction 34.9%
- Trend diagnostics: 19/29 NOAA sites significant (MK), 1/100 ADCIRC (chance level)
- Detrending robustness: fusion gains survive (34.7% vs 34.9% RL RMSE reduction)
- Package generalization: source_params abstraction across all core functions
- Package data (gulf_data.rda), README.md, vignette
- Regression check: 18/18 pass, zero numerical change from refactoring
- 33 tests passing

### Key results
- Cross-source correlations: mu 0.992, xi 0.952, log_sigma 0.521
- Fusion reduces return level SE by 35% overall vs NOAA-only
- LOO-CV: 68% reduction in mu prediction error, 35% in 100-yr RL RMSE
- Robust to detrending for sea level rise

### Saved model artifacts
- `data-raw/model_6dim_best.rds` — best joint model (NLL = -103.61)
- `data-raw/model_noaa_only.rds` / `model_adcirc_only.rds` — baselines
- `data-raw/loo_cv_results.rds` — LOO-CV results
- `data-raw/detrended_results.rds` — full detrended pipeline results
- `data-raw/trend_diagnostics.rds` — trend test results at 29 NOAA sites
- `data-raw/predictions_grid.rds` / `return_levels_100yr.rds` — predictions

## Next steps
1. Taper sensitivity (lambda in {75, 150, 300})
2. Manuscript preparation
