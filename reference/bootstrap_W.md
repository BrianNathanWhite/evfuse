# Estimate W via nonparametric block bootstrap (stationary)

Implements the bootstrap procedure from Section 2.4 of Russell et al.
(2020). Resamples years (preserving spatial dependence) and re-fits
stationary GEV at each site, then estimates the covariance of the
stage-1 estimators. For a nonstationary alternative that accounts for
NOAA sea level trends, see
[`bootstrap_W_detrended`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md).

## Usage

``` r
bootstrap_W(dat, B = 500, log_scale = TRUE, seed = NULL)
```

## Arguments

- dat:

  An `evfuse_data` object.

- B:

  Number of bootstrap replications (default 500).

- log_scale:

  Logical. If TRUE, work with log(sigma) (default TRUE).

- seed:

  Random seed for reproducibility.

## Value

A list with components:

- W_bs:

  The raw bootstrap covariance matrix (Lp x Lp).

- Gamma:

  The bootstrap matrix of MLEs (B x Lp).

- n_failures:

  Number of (site, bootstrap) combinations that failed.

## Note

A small fraction of bootstrap resamples may produce degenerate data that
causes warnings or convergence failures in
[`extRemes::fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html). Failed
fits are recorded in `n_failures` and excluded via pairwise-complete
covariance estimation (failure rates are typically well below 1
percent).

## See also

[`bootstrap_W_detrended`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md)
for the nonstationary version.
