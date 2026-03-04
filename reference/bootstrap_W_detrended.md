# Bootstrap W with NOAA detrending

Like
[`bootstrap_W`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W.md)
but NOAA sites are fitted with nonstationary GEV (linear mu trend),
extracting (mu0, log_sigma, xi) per replicate. ADCIRC sites use the
standard stationary GEV. This is the primary bootstrap approach used in
the manuscript.

## Usage

``` r
bootstrap_W_detrended(dat, df, B = 500, ref_year = 2000, seed = NULL)
```

## Arguments

- dat:

  An `evfuse_data` object.

- df:

  The raw data frame with year column.

- B:

  Number of bootstrap replications.

- ref_year:

  Reference year for centering.

- seed:

  Random seed.

## Value

A list with W_bs, Gamma, and n_failures.

## Note

Some bootstrap resamples will produce degenerate data that fails GEV
fitting (convergence warnings from
[`extRemes::fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html)). These
are recorded in `n_failures` and excluded via pairwise-complete
covariance estimation. Failure rates are typically well under 1 percent.

## See also

[`bootstrap_W`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W.md)
for the stationary version.
