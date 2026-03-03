# Gradient of return level w.r.t. (mu, log_sigma, xi)

For the delta method. Note: parameterized in terms of log(sigma), not
sigma, to match the kriging output.

## Usage

``` r
rl_gradient(mu, log_sigma, xi, r)
```

## Arguments

- mu:

  Location parameter.

- log_sigma:

  Log of scale parameter.

- xi:

  Shape parameter.

- r:

  Return period.

## Value

Numeric vector of length 3.
