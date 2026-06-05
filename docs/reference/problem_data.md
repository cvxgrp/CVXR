# Get Problem Data for a Solver

Returns the problem data that would be passed to a specific solver,
along with the reduction chain and inverse data for solution retrieval.

## Usage

``` r
problem_data(
  x,
  solver = NULL,
  gp = FALSE,
  enforce_dpp = FALSE,
  ignore_dpp = FALSE,
  ...
)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- solver:

  Character string naming solver, or `NULL` for automatic selection.

- gp:

  Logical; if `TRUE`, parse the problem as a geometric program.

- enforce_dpp:

  Logical; if `TRUE`, raise an error when a parametrized problem is not
  DPP instead of compiling it as non-DPP.

- ignore_dpp:

  Logical; if `TRUE`, treat a DPP problem as non-DPP (skip the DPP fast
  path).

- ...:

  Additional solver options.

## Value

A list with components `data`, `chain`, and `inverse_data`.
