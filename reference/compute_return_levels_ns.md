# Compute nonstationary return levels with covariate

Takes kriging predictions of 4 GEV parameters (mu0, mu1, log_sigma, xi)
and a future covariate value, and computes r-year return levels with
confidence intervals via the delta method and/or simulation.

## Usage

``` r
compute_return_levels_ns(
  predictions,
  covariate_value,
  r = 100,
  alpha = 0.05,
  method = "both",
  n_sim = 2500,
  seed = NULL
)
```

## Arguments

- predictions:

  An `evfuse_predictions` object from
  [`predict_krig`](https://briannathanwhite.github.io/evfuse/reference/predict_krig.md)
  with 4-column `noaa_mean` and 4x4 covariance matrices.

- covariate_value:

  Scalar or vector (length n_new) of centered covariate values at
  prediction sites.

- r:

  Return period in years (default 100).

- alpha:

  Confidence level for CIs (default 0.05 for 95% CIs).

- method:

  One of `"delta"`, `"simulation"`, or `"both"` (default).

- n_sim:

  Number of simulations (default 2500).

- seed:

  Random seed for simulation method.

## Value

A data frame with columns: lon, lat, return_level, se_delta,
ci_lower_delta, ci_upper_delta, se_sim, ci_lower_sim, ci_upper_sim.

## Details

The effective location is `mu = mu0 + mu1 * covariate_value`, where
`covariate_value` should be centered relative to the same reference used
in
[`fit_gev_ns`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_ns.md).
