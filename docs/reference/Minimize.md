# Create a Minimization Objective

Specifies that the objective expression should be minimized. The
expression must be convex and scalar for a DCP-compliant problem.

## Usage

``` r
Minimize(expr)
```

## Arguments

- expr:

  A CVXR expression or numeric value to minimize.

## Value

A `Minimize` object.

## See also

[`Maximize`](https://www.cvxgrp.org/CVXR/reference/Maximize.md),
[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)

## Examples

``` r
x <- Variable()
obj <- Minimize(x^2 + 1)
```
