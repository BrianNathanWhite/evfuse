# Predict at new locations via universal kriging

Given the fitted stage 2 model, predict the full parameter vector and
its covariance at new spatial locations. Then extract the target
components (default: first source's params, e.g. 1-3 for NOAA) and their
marginal covariance.

## Usage

``` r
predict_krig(model, new_sites, predict_params = NULL)
```

## Arguments

- model:

  A fitted `evfuse_model` from `fit_spatial_model`.

- new_sites:

  Data frame with lon and lat columns for prediction locations.

- predict_params:

  Integer vector of parameter indices to extract for predictions
  (default: first source's params, i.e. `1:3`).

## Value

A list with components:

- pred_mean:

  Matrix (n_new x p) of predicted parameter means.

- pred_cov:

  List of n_new pxp predictive covariance matrices.

- noaa_mean:

  Matrix (n_new x length(predict_params)) of predicted target GEV
  parameters (mu, log_sigma, xi).

- noaa_cov:

  List of n_new marginal covariances for target params.
