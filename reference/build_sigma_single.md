# Marginal covariance at a single site (no spatial separation)

At distance 0, Omega_i = 1, so Sigma = A A^T.

## Usage

``` r
build_sigma_single(A, rho, p = 6)
```

## Arguments

- A:

  Lower triangular matrix.

- rho:

  Range parameters (unused at distance 0, included for API consistency).

- p:

  Dimension.

## Value

p x p covariance matrix.
