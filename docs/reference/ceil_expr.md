# Elementwise Ceiling

Returns the ceiling (smallest integer \>= x) of each element. This atom
is quasiconvex and quasiconcave but NOT convex or concave, so it can
only be used in DQCP problems (solved with `qcp = TRUE`).

## Usage

``` r
ceil_expr(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

A `Ceil` expression.

## See also

[`floor_expr`](https://www.cvxgrp.org/CVXR/reference/floor_expr.md)
