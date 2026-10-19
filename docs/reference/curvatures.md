# All Curvatures of an Expression

Lists every curvature that holds for an expression, DCP, DGP and DQCP
alike, as CVXPY's `Expression.curvatures`.
[`curvature`](https://www.cvxgrp.org/CVXR/reference/curvature.md)
reports only the strongest DCP curvature.

## Usage

``` r
curvatures(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

A character vector, in the order `"CONSTANT"`, `"AFFINE"`, `"CONVEX"`,
`"CONCAVE"`, `"LOG-LOG CONSTANT"`, `"LOG-LOG AFFINE"`,
`"LOG-LOG CONVEX"`, `"LOG-LOG CONCAVE"`, `"QUASILINEAR"`,
`"QUASICONVEX"`, `"QUASICONCAVE"`, holding those that apply; `"UNKNOWN"`
alone if none does.

## See also

[`curvature()`](https://www.cvxgrp.org/CVXR/reference/curvature.md),
[`log_log_curvature()`](https://www.cvxgrp.org/CVXR/reference/log_log_curvature.md)
