# Embed W_tap into full parameter space

The bootstrap produces a covariance matrix in observed-parameter space
(L sites x p_bs params per site). This function embeds it into the full
joint parameter space (L sites x p_full params) required by the stage 2
model. Each source's observed params are mapped to their corresponding
positions in the joint model via `source_params`.

## Usage

``` r
embed_W(W_tap, dat, source_params = NULL)
```

## Arguments

- W_tap:

  Tapered bootstrap covariance (L*p_bs x L*p_bs).

- dat:

  Data object with site information and `source_params`.

- source_params:

  Named list mapping source labels to full-model param indices. Default:
  `list(NOAA = 1:3, ADCIRC = 4:6)`.

## Value

Embedded covariance matrix (L*p_full x L*p_full).

## Details

When sources observe different numbers of parameters (e.g., NOAA has 4
and ADCIRC has 3 in the nonstationary case), the bootstrap uses a union
parameter ordering with `NA` entries for unobserved positions. This
function detects unobserved positions from `NA` on the diagonal of
`W_tap` and maps only the observed entries.
