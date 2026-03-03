# GEV trend test via likelihood ratio

Fits stationary GEV(mu, sigma, xi) and nonstationary GEV(mu0 +
mu1\*year, sigma, xi) via
[`extRemes::fevd`](https://rdrr.io/pkg/extRemes/man/fevd.html), then
compares with a likelihood ratio test (chi-squared, df=1).

## Usage

``` r
gev_trend_test(annual_maxima, years)
```

## Arguments

- annual_maxima:

  Numeric vector of annual maxima.

- years:

  Numeric vector of corresponding years (same length).

## Value

A list with components:

- mu0:

  Intercept of nonstationary location.

- mu1:

  Slope of nonstationary location (units per year).

- lrt_stat:

  Likelihood ratio test statistic.

- lrt_pvalue:

  P-value from chi-squared(df=1).

- aic_stat:

  AIC of stationary model.

- aic_nonstat:

  AIC of nonstationary model.

- nllh_stat:

  Negative log-likelihood of stationary model.

- nllh_nonstat:

  Negative log-likelihood of nonstationary model.
