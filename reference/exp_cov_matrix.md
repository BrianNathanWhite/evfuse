# Exponential covariance matrix

Computes the correlation matrix Omega_i where Omega_i(j,k) = exp(-D(j,k)
/ rho) for range parameter rho \> 0.

## Usage

``` r
exp_cov_matrix(D, rho)
```

## Arguments

- D:

  Distance matrix.

- rho:

  Range parameter (positive).

## Value

Correlation matrix of same dimension as D.
