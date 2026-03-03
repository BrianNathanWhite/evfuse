# Compute return levels with confidence intervals

Takes kriging predictions (NOAA GEV parameters and their covariance) and
computes r-year return levels with 95% confidence intervals via the
delta method and/or simulation from the predictive distribution.

## Usage

``` r
compute_return_levels(
  predictions,
  r = 100,
  alpha = 0.05,
  method = "both",
  n_sim = 2500,
  seed = NULL
)
```

## Arguments

- predictions:

  An `evfuse_predictions` object from `predict_krig`.

- r:

  Return period in years (default 100).

- alpha:

  Confidence level for CIs (default 0.05 for 95% CIs).

- method:

  One of "delta", "simulation", or "both" (default "both").

- n_sim:

  Number of simulations for the simulation method (default 2500).

- seed:

  Random seed for simulation method.

## Value

A data frame with columns:

- lon, lat:

  Coordinates of prediction sites.

- return_level:

  Point estimate of the r-year return level.

- se_delta:

  Standard error from the delta method (if computed).

- ci_lower_delta, ci_upper_delta:

  Delta method CI bounds.

- se_sim:

  Standard error from simulation (if computed).

- ci_lower_sim, ci_upper_sim:

  Simulation-based CI bounds.
