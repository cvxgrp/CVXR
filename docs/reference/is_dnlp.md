# Check if an Expression or Problem is DNLP-Compliant

Tests whether the object follows the Disciplined Nonlinear Programming
(DNLP) rules: smooth representable, i.e. linearizable-convex or
linearizable-concave.

## Usage

``` r
is_dnlp(x, ...)
```

## Arguments

- x:

  An expression, objective, constraint, or
  [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md).

- ...:

  Not used.

## Value

Logical scalar.
