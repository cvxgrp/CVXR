# Create an Optimization Problem

Constructs a convex optimization problem from an objective and a list of
constraints. Use
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md) to solve the
problem.

## Usage

``` r
Problem(objective, constraints = list())
```

## Arguments

- objective:

  A [`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md) or
  [`Maximize`](https://www.cvxgrp.org/CVXR/reference/Maximize.md)
  object.

- constraints:

  A list of `Constraint` objects (e.g., created by `==`, `<=`, `>=`
  operators on expressions). Defaults to an empty list (unconstrained).

## Value

A `Problem` object.

## Known limitations

- Problems must contain at least one
  [`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md).
  Zero-variable problems (e.g., minimizing a constant) will cause an
  internal error in the reduction pipeline.

## Examples

``` r
x <- Variable(2)
prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
```
