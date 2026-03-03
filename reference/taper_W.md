# Apply covariance tapering to W

Computes W_tap = W_bs \* T_tap (Hadamard product) where T_tap is built
from the Wendland 2 function with range lambda.

## Usage

``` r
taper_W(W_bs, D, lambda, p = 3)
```

## Arguments

- W_bs:

  Raw bootstrap covariance matrix (Lp x Lp).

- D:

  Distance matrix (L x L).

- lambda:

  Taper range in km.

- p:

  Number of GEV parameters per site (default 3).

## Value

Tapered covariance matrix W_tap (Lp x Lp).
