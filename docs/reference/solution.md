# Get the Raw Solution Object

Returns the raw `Solution` object from the most recent solve, containing
primal and dual variable values, status, and solver attributes.

## Usage

``` r
solution(x)
```

## Arguments

- x:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

## Value

A `Solution` object, or `NULL` if the problem has not been solved.
