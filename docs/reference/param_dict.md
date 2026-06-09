# Get all Parameters of a Problem as a Named List

Mirrors CVXPY's `Problem.param_dict` property
(`cvxpy/problems/problem.py:260-264`): returns a named list keyed by
each parameter's name, where the value is the
[`Parameter`](https://www.cvxgrp.org/CVXR/reference/Parameter.md) object
itself.

## Usage

``` r
param_dict(x, ...)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- ...:

  Not used.

## Value

Named list of
[`Parameter`](https://www.cvxgrp.org/CVXR/reference/Parameter.md)
objects, keyed by name.
