# Build observation structure for partial observations

Maps each site to the indices of its observed components in the full
Lp-dimensional parameter vector.

## Usage

``` r
build_observation_structure(stage1, dat)
```

## Arguments

- stage1:

  Stage 1 fits.

- dat:

  Data object.

## Value

List with obs_idx (indices into the full vector) and theta_obs (the
observed parameter values).
