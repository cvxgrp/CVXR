# Elementwise Floor

Returns the floor (largest integer \<= x) of each element. This atom is
quasiconvex and quasiconcave but NOT convex or concave, so it can only
be used in DQCP problems (solved with `qcp = TRUE`). As in CVXPY, its
value rounds `x` to 4 decimal places first, so `floor(-3e-10)` of an
expression is 0, not -1.

## Usage

``` r
floor_expr(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

A `Floor` expression.

## See also

[`ceil_expr`](https://www.cvxgrp.org/CVXR/reference/ceil_expr.md)
