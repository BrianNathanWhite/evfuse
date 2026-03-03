# Package index

## Package

- [`evfuse`](https://briannathanwhite.github.io/evfuse/reference/evfuse-package.md)
  [`evfuse-package`](https://briannathanwhite.github.io/evfuse/reference/evfuse-package.md)
  : evfuse: Spatial Data Fusion for Extreme Value Analysis

## Data

Load data and compute distance matrices

- [`coast_data`](https://briannathanwhite.github.io/evfuse/reference/coast_data.md)
  : U.S. coastal sea level data
- [`prediction_grid`](https://briannathanwhite.github.io/evfuse/reference/prediction_grid.md)
  : Coastal prediction grid
- [`load_data()`](https://briannathanwhite.github.io/evfuse/reference/load_data.md)
  : Load and validate sea level data
- [`compute_distances()`](https://briannathanwhite.github.io/evfuse/reference/compute_distances.md)
  : Compute pairwise great-circle distances between sites
- [`compute_cross_distances()`](https://briannathanwhite.github.io/evfuse/reference/compute_cross_distances.md)
  : Cross-distance matrix between two sets of sites

## Stage 1: GEV Fitting

Pointwise generalized extreme value fits

- [`fit_gev_all()`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_all.md)
  : Fit GEV at all sites (Stage 1, stationary)
- [`fit_gev_detrended()`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md)
  : Fit GEV at all sites with NOAA detrending (Stage 1, nonstationary)
- [`fit_gev_ns()`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_ns.md)
  : Fit nonstationary GEV with covariate (Stage 1)
- [`gev_return_level()`](https://briannathanwhite.github.io/evfuse/reference/gev_return_level.md)
  : GEV quantile function (return level)
- [`gev_exceedance_prob()`](https://briannathanwhite.github.io/evfuse/reference/gev_exceedance_prob.md)
  : GEV exceedance probability

## Bootstrap

Measurement uncertainty via block bootstrap

- [`bootstrap_W()`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W.md)
  : Estimate W via nonparametric block bootstrap (stationary)
- [`bootstrap_W_detrended()`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md)
  : Bootstrap W with NOAA detrending
- [`bootstrap_W_ns()`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_ns.md)
  : Bootstrap W with nonstationary covariate
- [`taper_W()`](https://briannathanwhite.github.io/evfuse/reference/taper_W.md)
  : Apply covariance tapering to W
- [`subset_W_bs()`](https://briannathanwhite.github.io/evfuse/reference/subset_W_bs.md)
  : Subset bootstrap covariance to a single data source
- [`embed_W()`](https://briannathanwhite.github.io/evfuse/reference/embed_W.md)
  : Embed W_tap into full parameter space

## Covariance

Spatial covariance and tapering functions

- [`exp_cov_matrix()`](https://briannathanwhite.github.io/evfuse/reference/exp_cov_matrix.md)
  : Exponential covariance matrix
- [`wendland2()`](https://briannathanwhite.github.io/evfuse/reference/wendland2.md)
  : Wendland 2 compactly supported covariance function
- [`wendland_taper_matrix()`](https://briannathanwhite.github.io/evfuse/reference/wendland_taper_matrix.md)
  : Wendland 2 taper correlation matrix
- [`build_sigma()`](https://briannathanwhite.github.io/evfuse/reference/build_sigma.md)
  : Build the coregionalization covariance matrix
- [`build_taper()`](https://briannathanwhite.github.io/evfuse/reference/build_taper.md)
  : Build the taper matrix for W

## Stage 2: Spatial Model

Joint coregionalization model fitting

- [`fit_spatial_model()`](https://briannathanwhite.github.io/evfuse/reference/fit_spatial_model.md)
  : Fit Stage 2 spatial model via MLE
- [`fit_naive_model()`](https://briannathanwhite.github.io/evfuse/reference/fit_naive_model.md)
  : Fit a naive (single-source) 3-dim spatial model
- [`verify_gradient()`](https://briannathanwhite.github.io/evfuse/reference/verify_gradient.md)
  : Verify analytic gradient against numerical finite differences

## Kriging & Prediction

Spatial prediction and return level estimation

- [`predict_krig()`](https://briannathanwhite.github.io/evfuse/reference/predict_krig.md)
  : Predict at new locations via universal kriging
- [`predict_krig_naive()`](https://briannathanwhite.github.io/evfuse/reference/predict_krig_naive.md)
  : Predict at new locations from a naive model
- [`compute_return_levels()`](https://briannathanwhite.github.io/evfuse/reference/compute_return_levels.md)
  : Compute return levels with confidence intervals
- [`compute_return_levels_ns()`](https://briannathanwhite.github.io/evfuse/reference/compute_return_levels_ns.md)
  : Compute nonstationary return levels with covariate

## Validation

Leave-one-out cross-validation

- [`loo_cv()`](https://briannathanwhite.github.io/evfuse/reference/loo_cv.md)
  : Leave-one-out cross-validation via Rasmussen & Williams (2006)
  shortcut
- [`loo_summary()`](https://briannathanwhite.github.io/evfuse/reference/loo_summary.md)
  : Summarize LOO-CV results

## Trend Diagnostics

Temporal trend testing

- [`mann_kendall_test()`](https://briannathanwhite.github.io/evfuse/reference/mann_kendall_test.md)
  : Mann-Kendall trend test with Sen's slope
- [`gev_trend_test()`](https://briannathanwhite.github.io/evfuse/reference/gev_trend_test.md)
  : GEV trend test via likelihood ratio
