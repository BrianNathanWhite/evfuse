# Fit Stage 2 spatial model via MLE

Estimates beta (6x1), A (6x6 lower triangular), and rho (6x1) by
maximizing the Gaussian likelihood from Eq. (7) of Russell et al.
(2020), extended to handle partial observations. NOAA sites observe
components 1-3 and ADCIRC sites observe components 4-6 of the
6-dimensional parameter vector.

## Usage

``` r
fit_spatial_model(
  stage1,
  dat,
  W_tap,
  D,
  start = NULL,
  method = "L-BFGS-B",
  control = list(maxit = 1000, trace = 1),
  check_gradient = FALSE
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

- W_tap:

  Tapered covariance matrix from `taper_W`.

- D:

  Distance matrix from `compute_distances`.

- start:

  Optional named list of starting values for beta, A, rho.

- method:

  Optimization method for `optim` (default "L-BFGS-B").

- control:

  Control list passed to `optim`.

- check_gradient:

  If TRUE, verify analytic gradient against numerical finite differences
  at the starting values before optimizing.

## Value

A list with components:

- beta:

  Estimated mean vector (6x1).

- A:

  Estimated lower triangular matrix (6x6).

- rho:

  Estimated range parameters (6x1).

- Sigma:

  Estimated covariance matrix at observed sites.

- optim_result:

  Raw output from optim.

- grad_check:

  If `check_gradient = TRUE`, the output of `verify_gradient`.

## Details

Uses an analytic gradient derived from the coregionalization structure
(see Sections 2, 5, 6 of the gradient notes). The key identity is that
dSigma/dA_ab and dSigma/drho_k are rank-1 Kronecker products, enabling
O(n_obs^2) gradient computation per parameter instead of O(n_obs^3).
