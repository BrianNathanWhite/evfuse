# Bootstrap W with nonstationary covariate

Like
[`bootstrap_W_detrended`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md)
but uses an arbitrary covariate (e.g., SST) instead of a linear time
trend, and retains the covariate sensitivity mu1 in the bootstrap
parameter vector. NOAA sites are fitted with nonstationary GEV; ADCIRC
sites use stationary GEV.

## Usage

``` r
bootstrap_W_ns(dat, df, covariate, B = 500, ref_value = NULL, seed = NULL)
```

## Arguments

- dat:

  An `evfuse_data` object.

- df:

  The raw data frame with covariate column.

- covariate:

  String column name in `df`.

- B:

  Number of bootstrap replications.

- ref_value:

  Centering value (default: mean across NOAA site-years).

- seed:

  Random seed.

## Value

A list with `W_bs` (4L x 4L), `Gamma` (B x 4L), and `n_failures`.

## Details

The bootstrap Gamma matrix has `L * 4` columns in parameter-major
ordering: (mu0, mu1, log_sigma, xi) at all sites. ADCIRC mu1 entries are
`NA`; the resulting `W_bs` has `NA` in those positions. This is by
design:
[`embed_W`](https://briannathanwhite.github.io/evfuse/reference/embed_W.md)
skips unobserved entries.

## Note

A small fraction of bootstrap resamples may produce degenerate data that
causes warnings or convergence failures in
[`extRemes::fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html). Failed
fits are recorded in `n_failures` and excluded via pairwise-complete
covariance estimation (failure rates are typically well below 1
percent).

## See also

[`bootstrap_W_detrended`](https://briannathanwhite.github.io/evfuse/reference/bootstrap_W_detrended.md)
for the time-trend version.
