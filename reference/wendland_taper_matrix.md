# Wendland 2 taper correlation matrix

Constructs a sparse taper matrix based on the Wendland C4 function with
compact support at distance lambda.

## Usage

``` r
wendland_taper_matrix(D, lambda)
```

## Arguments

- D:

  Distance matrix.

- lambda:

  Taper range (same units as D). Correlations are zero beyond this
  distance.

## Value

Sparse taper matrix (class "dgCMatrix" from Matrix package).
