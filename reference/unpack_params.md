# Unpack parameter vector into beta, A, rho

The last p entries are log(rho); this function exponentiates to return
rho in the natural (positive) scale.

## Usage

``` r
unpack_params(par, p = 6)
```

## Arguments

- par:

  Numeric vector from optim.

- p:

  Dimension (default 6).

## Value

Named list with beta, A, rho (positive scale).
