# Summarize LOO-CV results

Computes RMSE, MAD, and return level metrics from LOO predictions.

## Usage

``` r
loo_summary(loo, r = 100)
```

## Arguments

- loo:

  An `evfuse_loo` object.

- r:

  Return period for return level comparison (default 100).

## Value

A list with parameter-level and return-level accuracy metrics.
