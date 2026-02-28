# Get Problem Data for a Solver

Returns the problem data that would be passed to a specific solver,
along with the reduction chain and inverse data for solution retrieval.

## Usage

``` r
problem_data(x, solver = NULL, ...)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- solver:

  Character string naming solver, or `NULL` for automatic selection.

- ...:

  Additional arguments.

## Value

A list with components `data`, `chain`, and `inverse_data`.
