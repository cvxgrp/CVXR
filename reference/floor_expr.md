# Elementwise Floor

Returns the floor (largest integer \<= x) of each element. This atom is
quasiconvex and quasiconcave but NOT convex or concave, so it can only
be used in DQCP problems (solved with `qcp = TRUE`).

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
