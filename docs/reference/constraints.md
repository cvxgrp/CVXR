# Get Problem Constraints (read-only)

Returns a copy of the problem's constraint list.

## Usage

``` r
constraints(x)
```

## Arguments

- x:

  A [Problem](https://www.cvxgrp.org/CVXR/reference/Problem.md) object.

## Value

A list of constraint objects.

## Details

Problem objects are **immutable**: constraints cannot be modified after
construction. To change constraints, create a new
[`Problem()`](https://www.cvxgrp.org/CVXR/reference/Problem.md). This
matches CVXPY's design where problems are immutable except through
[Parameter](https://www.cvxgrp.org/CVXR/reference/Parameter.md) value
changes.

## See also

[`Problem()`](https://www.cvxgrp.org/CVXR/reference/Problem.md),
[`objective()`](https://www.cvxgrp.org/CVXR/reference/objective.md)

## Examples

``` r
x <- Variable(2)
prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
length(constraints(prob))  # 1
#> [1] 1
```
