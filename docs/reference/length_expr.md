# Length of a Vector (Last Nonzero Index)

Returns the index of the last nonzero element of a vector (1-based).
This atom is quasiconvex but NOT convex, so it can only be used in DQCP
problems (solved with `qcp = TRUE`).

## Usage

``` r
length_expr(x)
```

## Arguments

- x:

  A CVXR vector expression.

## Value

A `Length` expression (scalar).

## See also

[`ceil_expr`](https://www.cvxgrp.org/CVXR/reference/ceil_expr.md),
[`floor_expr`](https://www.cvxgrp.org/CVXR/reference/floor_expr.md)
