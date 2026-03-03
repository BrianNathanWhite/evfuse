# Fit GEV at all sites with NOAA detrending (Stage 1, nonstationary)

NOAA sites get nonstationary GEV: mu(t) = mu0 + mu1\*(t - ref_year),
sigma and xi stationary. The detrended triplet (mu0, log_sigma, xi) is
returned. ADCIRC sites get the usual stationary fit. This is the primary
Stage 1 approach used in the manuscript.

## Usage

``` r
fit_gev_detrended(dat, df, ref_year = 2000)
```

## Arguments

- dat:

  An `evfuse_data` object.

- df:

  The raw data frame with year column.

- ref_year:

  Reference year for centering (default 2000).

## Value

An `evfuse_stage1` object with detrended NOAA parameters.

## See also

[`fit_gev_all`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_all.md)
for the stationary version.
