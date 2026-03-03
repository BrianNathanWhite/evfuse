# Verify analytic gradient against numerical finite differences

Computes the analytic gradient and a central finite-difference
approximation at the same parameter vector, then reports the relative
error for each component.

## Usage

``` r
verify_gradient(par, fn, gr, eps = 1e-05)
```

## Arguments

- par:

  Parameter vector to check at.

- fn:

  Objective function (scalar-valued).

- gr:

  Gradient function (returns vector same length as par).

- eps:

  Finite difference step size (default 1e-5).

## Value

List with components:

- analytic:

  Analytic gradient vector.

- numeric:

  Numerical gradient vector.

- rel_err:

  Relative error for each component.

- max_rel_err:

  Maximum relative error across all components.
