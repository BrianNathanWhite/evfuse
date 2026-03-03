# Leave-one-out cross-validation via Rasmussen & Williams (2006) shortcut

Computes closed-form LOO-CV predictions at NOAA sites using the block
version of Rasmussen & Williams (2006) eq. 5.12. For site k with
observation block B_k, the LOO predictive mean and covariance are
computed from the precision matrix without refitting.

## Usage

``` r
loo_cv(model, ...)

# S3 method for class 'evfuse_model'
loo_cv(model, loo_source = NULL, ...)

# S3 method for class 'evfuse_naive_model'
loo_cv(model, ...)
```

## Arguments

- model:

  A fitted model (evfuse_model or evfuse_naive_model).

- ...:

  Additional arguments (currently unused).

- loo_source:

  Which data source to evaluate LOO-CV at (default: first source in
  `source_params`, typically `"NOAA"`).

## Value

An `evfuse_loo` object with components:

- loo_mean:

  Matrix (n_sites x 3) of LOO predictive means.

- loo_cov:

  List of 3x3 LOO predictive covariance matrices.

- observed:

  Matrix (n_sites x 3) of observed Stage 1 MLEs.

- sites:

  Data frame of site metadata.

- loo_lpd:

  Vector of per-site log predictive densities.
