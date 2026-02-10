# evfuse

Spatial data fusion for extreme value analysis of coastal sea levels along the US East and Gulf Coasts. Implements the Russell et al. (2019) two-stage frequentist framework extended to 6-dimensional coregionalization with partial observations, fusing NOAA tidal gauge data (29 sites) with ADCIRC storm surge simulations (100 sites).

## Key conventions
- Parameter-major ordering throughout (parameter index varies fastest)
- Log-rho parameterization for optimization (ensures rho > 0)
- Lambda = 300 km Wendland C2 taper for bootstrap W
- B = 500 bootstrap replicates
- 20 multi-start random restarts from log-uniform[50, 5000] km
- Prediction grid: filtered coastal points at 15 km spacing, lon -98 to -66, lat 24 to 46
- Seeds: bootstrap (42), retry loop (101-105), multi-start (2026), return level simulation (123)

## Project structure
- R/data.R — data loading and distance matrix
- R/gev.R — Stage 1 pointwise GEV fitting via extRemes::fevd
- R/covariance.R — exponential covariance, Wendland taper, Kronecker construction
- R/bootstrap.R — block bootstrap for measurement error covariance W
- R/spatial_model.R — Stage 2 MLE with analytic gradient (6-dim joint model)
- R/naive_model.R — 3-dim baseline models (NOAA-only, ADCIRC-only)
- R/kriging.R — GP kriging with partial observations
- R/return_levels.R — 100-year return levels via delta method and simulation
- scripts/run_pipeline.R — reproducible full pipeline (Stage 1 -> bootstrap -> Stage 2 -> krige -> return levels)
- scripts/run_naive.R — fit and evaluate naive baseline models
- scripts/multi_start.R — multi-start optimization for Stage 2
- scripts/plot_results.R — main visualization (Stage 1 maps, kriged params, return levels)
- scripts/plot_comparison.R — comparison plots and tables (joint vs NOAA-only vs ADCIRC-only)
- scripts/print_tables.R — formatted console output of comparison tables
- scripts/save_tables.R — save comparison tables as text and RDS
- data-raw/ — source data (.csv), fitted models (.rds), prediction grid, comparison tables
- figures/ — all output plots
- tests/testthat/test-core.R — 33 unit tests
- evfuse_notes.tex — running working notes with analytical descriptions of all results

## Data
- Source: `data-raw/combined_max.csv`
- 129 sites: 29 NOAA, 100 ADCIRC
- 34-43 years of annual maxima per site

## Development practices
- Use `devtools::load_all()` and `devtools::test()` to iterate
- Run `test()` after every change
- Don't install packages without asking first
- If a fix involves changing mathematical structure (not just R syntax), explain what's changing and why before implementing

## Implementation notes
- `extRemes::fevd` returns observed information matrix directly (positive definite), not Hessian of negative log-likelihood. Use `solve(fit$results$hessian)` directly for vcov.
- Delta method for log_sigma transformation: J = diag(1, 1/sigma, 1)
- Rho optimized in log space: `pack_params` stores `log(rho)`, `unpack_params` exponentiates. Gradient uses chain rule.
- Analytic gradient contracts Q = V_inv - alpha*alpha^T with Omega_i into 6x6 H matrices. See `data-raw/evfuse_gradient.pdf`.
- Sigma_obs built directly in 387x387 observed space (never form full 774x774).
- nll/nll_grad share a cache keyed on parameter vector to avoid redundant Cholesky.

## Current state (Feb 2026)
- Stage 1: All 129/129 sites converge, all vcov matrices PSD
- Stage 2: Joint 6-dim model fitted (NLL = -103.61, 33 params, lambda=300, B=500)
- Naive baselines: NOAA-only (NLL = +13.67, 12 params) and ADCIRC-only (NLL = -75.49, 12 params)
- Kriging predictions computed on coastal prediction grid for all three models
- 100-year return levels with delta-method SE and simulation CIs for all three models
- Comparison tables saved to data-raw/comparison_tables.rds and data-raw/comparison_tables.txt
- All diagnostic and comparison figures generated in figures/
- 33 tests passing

### Saved model artifacts
- `data-raw/model_6dim_best.rds` — best joint model from multi-start (NLL = -103.61)
- `data-raw/model_noaa_only.rds` — NOAA-only 3-dim baseline
- `data-raw/model_adcirc_only.rds` — ADCIRC-only 3-dim baseline
- `data-raw/predictions_grid.rds` — joint kriging predictions on grid
- `data-raw/return_levels_100yr.rds` — joint 100-yr return levels
- `data-raw/rl_noaa_only.rds`, `data-raw/rl_adcirc_only.rds` — baseline return levels

## Key results
- Cross-source correlations: mu 0.992, xi 0.952, log_sigma 0.521
- Fusion reduces return level SE by 35% overall vs NOAA-only
- Median SE ratio (NOAA-only / joint) = 1.59; 57% of coast with ratio > 1.5
- Joint model improves xi RMSE at NOAA sites (0.111 vs 0.124)

## Next steps
1. LOO-CV at 29 NOAA sites (closed-form Rasmussen & Williams shortcut) for joint vs NOAA-only
2. Trend diagnostics (Mann-Kendall + GEV with linear mu(t)) to justify stationarity assumption
3. Taper sensitivity (lambda in {75, 150, 300})
4. Manuscript for Environmetrics
