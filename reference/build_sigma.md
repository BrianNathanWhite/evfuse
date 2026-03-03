# Build the coregionalization covariance matrix

Constructs the covariance matrix for the spatial random effects
following Eq. (7) in Russell et al. (2020). For p parameters and L
sites, Sigma = (I_L x A) sum_i (e_i e_i^T x Omega_i(rho)) (I_L x A)^T.

## Usage

``` r
build_sigma(A, rho, D, p = 6)
```

## Arguments

- A:

  Lower triangular matrix (p x p).

- rho:

  Vector of p range parameters.

- D:

  Distance matrix (L x L).

- p:

  Number of parameters (default 6).

## Value

Covariance matrix of dimension (L*p x L*p).

## Details

In our setting, p = 6 (3 NOAA + 3 ADCIRC GEV parameters).
