# Subset bootstrap covariance to a single data source

Extracts the source-specific block from the full bootstrap covariance
matrix W_bs (L*3 x L*3 in parameter-major ordering).

## Usage

``` r
subset_W_bs(W_bs, dat, source)
```

## Arguments

- W_bs:

  Raw bootstrap covariance matrix (L*3 x L*3).

- dat:

  An `evfuse_data` object.

- source:

  Character: "NOAA" or "ADCIRC".

## Value

Subsetted covariance matrix (L_sub*3 x L_sub*3).
