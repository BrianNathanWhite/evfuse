# Mann-Kendall trend test with Sen's slope

Tests for a monotonic trend in a time series using Kendall's tau
statistic. Implements the test directly from rank correlations (no
external package required). Sen's slope is the median of all pairwise
slopes.

## Usage

``` r
mann_kendall_test(annual_maxima, years)
```

## Arguments

- annual_maxima:

  Numeric vector of annual maxima.

- years:

  Numeric vector of corresponding years (same length).

## Value

A list with components:

- tau:

  Kendall's tau statistic.

- S:

  Mann-Kendall S statistic.

- p_value:

  Two-sided p-value (normal approximation).

- sens_slope:

  Sen's slope estimate (units per year).

- sens_intercept:

  Sen's intercept estimate.
