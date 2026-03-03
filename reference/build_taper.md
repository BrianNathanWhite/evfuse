# Build the taper matrix for W

Constructs T_tap(lambda) = 1_p 1_p^T kron C_W2(lambda) as in Russell et
al. This applies the same Wendland taper across all p parameter blocks,
allowing cross-parameter dependence at nearby stations.

## Usage

``` r
build_taper(D, lambda, p = 6)
```

## Arguments

- D:

  Distance matrix (L x L).

- lambda:

  Taper range in km.

- p:

  Number of GEV parameters (default 6).

## Value

Taper matrix of dimension (L*p x L*p).
