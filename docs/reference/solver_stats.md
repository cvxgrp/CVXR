# Get Solver Statistics

Returns solver statistics from the most recent solve, including solve
time, setup time, and iteration count.

## Usage

``` r
solver_stats(x)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

## Value

A `SolverStats` object, or `NULL` if the problem has not been solved.
