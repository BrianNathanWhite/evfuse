# Fit nonstationary GEV with covariate (Stage 1)

NOAA sites get nonstationary GEV: mu(t) = mu0 + mu1 \* (cov(t) -
ref_value), sigma and xi stationary. All four parameters (mu0, mu1,
log_sigma, xi) are returned for spatial modeling. ADCIRC sites get the
usual stationary fit with three parameters (mu0, log_sigma, xi).

## Usage

``` r
fit_gev_ns(dat, df, covariate, ref_value = NULL)
```

## Arguments

- dat:

  An `evfuse_data` object.

- df:

  The raw data frame containing the covariate column.

- covariate:

  String column name in `df` (e.g., `"sst_regional_warm"`).

- ref_value:

  Centering value for the covariate. Default: mean across all NOAA
  site-years. Centering improves numerical stability and makes mu0
  interpretable as the location at the reference covariate value.

## Value

An `evfuse_stage1` object with 4-column `theta_hat`. NOAA rows contain
(mu0, mu1, log_sigma, xi); ADCIRC rows contain (mu0, log_sigma, xi, NA).
The `vcov_list` entries are 4x4 for NOAA and 3x3 for ADCIRC.

## See also

[`fit_gev_detrended`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md)
for the time-trend version that discards mu1 before spatial modeling.
