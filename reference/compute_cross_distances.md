# Cross-distance matrix between two sets of sites

Cross-distance matrix between two sets of sites

## Usage

``` r
compute_cross_distances(sites1, sites2)
```

## Arguments

- sites1:

  Data frame with lon, lat (n1 rows).

- sites2:

  Data frame with lon, lat (n2 rows).

## Value

Matrix of dimension (n1 x n2) of great-circle distances in km.
