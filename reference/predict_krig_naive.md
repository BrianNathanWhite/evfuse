# Predict at new locations from a naive model

Universal kriging for the 3-dim naive model. All observed sites
contribute all 3 parameters (no partial observations).

## Usage

``` r
predict_krig_naive(model, new_sites)
```

## Arguments

- model:

  An `evfuse_naive_model` from `fit_naive_model`.

- new_sites:

  Data frame with lon and lat columns.

## Value

An `evfuse_predictions` object compatible with `compute_return_levels`.
