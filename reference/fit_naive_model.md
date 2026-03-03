# Fit a naive (single-source) 3-dim spatial model

Fits a 3-dimensional GP with coregionalization to data from a single
source (NOAA or ADCIRC only). All sites observe all 3 parameters,
eliminating the partial-observation complexity of the joint model.

## Usage

``` r
fit_naive_model(
  stage1,
  dat,
  W_bs,
  D,
  source = c("NOAA", "ADCIRC"),
  lambda = 300,
  n_starts = 20,
  control = list(maxit = 2000, trace = 1)
)
```

## Arguments

- stage1:

  An `evfuse_stage1` object from
  [`fit_gev_all`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_all.md)
  or
  [`fit_gev_detrended`](https://briannathanwhite.github.io/evfuse/reference/fit_gev_detrended.md).

- dat:

  An `evfuse_data` object.

- W_bs:

  Raw bootstrap covariance matrix (L*3 x L*3).

- D:

  Distance matrix (L x L).

- source:

  Character: "NOAA" or "ADCIRC".

- lambda:

  Wendland taper range in km (default 300).

- n_starts:

  Number of multi-start attempts (default 20).

- control:

  Control list passed to `optim`.

## Value

An `evfuse_naive_model` object.
