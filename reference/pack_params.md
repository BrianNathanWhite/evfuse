# Pack parameters into a single vector for optim

Rho is stored as log(rho) so the optimizer works in unconstrained space.

## Usage

``` r
pack_params(beta, A, rho, p = 6)
```

## Arguments

- beta:

  Mean vector (p).

- A:

  Lower triangular matrix (p x p).

- rho:

  Range parameters (p), in natural (positive) scale.

- p:

  Dimension (default 6).

## Value

Numeric vector.
