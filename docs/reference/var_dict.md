# Get all Variables of a Problem as a Named List

Mirrors CVXPY's `Problem.var_dict` property
(`cvxpy/problems/problem.py:267-271`): returns a named list keyed by
each variable's name, where the value is the
[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md) object
itself.

## Usage

``` r
var_dict(x, ...)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- ...:

  Not used.

## Value

Named list of
[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md) objects,
keyed by name.
