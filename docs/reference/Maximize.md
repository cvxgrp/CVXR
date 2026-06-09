# Create a Maximization Objective

Specifies that the objective expression should be maximized. The
expression must be concave and scalar for a DCP-compliant problem.

## Usage

``` r
Maximize(expr)
```

## Arguments

- expr:

  A CVXR expression or numeric value to maximize.

## Value

A `Maximize` object.

## See also

[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md),
[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)

## Examples

``` r
x <- Variable()
obj <- Maximize(-x^2 + 1)
```
