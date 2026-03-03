# Cross-covariance for a single new site

Builds the p x (L\*p) cross-covariance matrix between a single new site
and all L observed sites.

## Usage

``` r
build_cross_sigma_single(A, rho, d_vec, L, p = 6)
```

## Arguments

- A:

  Lower triangular matrix (p x p).

- rho:

  Range parameters (p).

- d_vec:

  Vector of distances from new site to each of L observed sites.

- L:

  Number of observed sites.

- p:

  Number of parameters (default 6).

## Value

Matrix of dimension (p x Lp).
