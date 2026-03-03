# Fit GEV at all sites (Stage 1, stationary)

Fits a stationary GEV distribution to annual block maxima at each site
using [`extRemes::fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html).
Returns pointwise MLEs and their asymptotic covariance matrices. For a
nonstationary alternative that accounts for sea level trends at NOAA
sites, see
[`fit_gev_detrended`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md).

## Usage

``` r
fit_gev_all(dat, log_scale = TRUE)
```

## Arguments

- dat:

  An `evfuse_data` object from `load_data`.

- log_scale:

  Logical. If TRUE (default), return log(sigma) instead of sigma. This
  is standard for the stage 2 model where we need unconstrained
  parameters.

## Value

A list with components:

- theta_hat:

  Matrix of dimension (n_sites x 3). Columns are mu, log_sigma (or
  sigma), xi. Rows ordered to match dat\$sites.

- vcov_list:

  List of n_sites 3x3 asymptotic covariance matrices for (mu,
  log_sigma, xi) at each site.

- fits:

  List of raw fevd fit objects for diagnostics.

- converged:

  Logical vector indicating convergence at each site.

## See also

[`fit_gev_detrended`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md)
for the nonstationary version.
