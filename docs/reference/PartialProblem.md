# Partial optimization of a Problem

A `PartialProblem` is an Expression that represents the optimal value of
an inner Problem as a function of the variables you choose NOT to
optimise over. Build one with
[`partial_optimize()`](https://www.cvxgrp.org/CVXR/reference/partial_optimize.md)
rather than constructing the class directly.

## Usage

``` r
PartialProblem(
  prob,
  opt_vars,
  dont_opt_vars,
  solver = NULL,
  solve_kwargs = list(),
  id = NULL
)
```
