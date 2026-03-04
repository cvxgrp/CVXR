# Get Problem Objective (read-only)

Returns the problem's objective.

## Usage

``` r
objective(x)
```

## Arguments

- x:

  A [Problem](https://www.cvxgrp.org/CVXR/reference/Problem.md) object.

## Value

A [Minimize](https://www.cvxgrp.org/CVXR/reference/Minimize.md) or
[Maximize](https://www.cvxgrp.org/CVXR/reference/Maximize.md) object.

## Details

Problem objects are **immutable**: the objective cannot be modified
after construction. To change the objective, create a new
[`Problem()`](https://www.cvxgrp.org/CVXR/reference/Problem.md).

## See also

[`Problem()`](https://www.cvxgrp.org/CVXR/reference/Problem.md),
[`constraints()`](https://www.cvxgrp.org/CVXR/reference/constraints.md)

## Examples

``` r
x <- Variable(2)
prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
objective(prob)
#> minimize SumEntries(var74, NULL, FALSE)
```
