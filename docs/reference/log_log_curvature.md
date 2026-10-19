# Log-Log Curvature of an Expression

The DGP analogue of
[`curvature`](https://www.cvxgrp.org/CVXR/reference/curvature.md), as
CVXPY's `Expression.log_log_curvature`.

## Usage

``` r
log_log_curvature(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

Character: the first of `"LOG-LOG CONSTANT"`, `"LOG-LOG AFFINE"`,
`"LOG-LOG CONVEX"`, `"LOG-LOG CONCAVE"` that holds, else `"UNKNOWN"`.

## See also

[`curvatures()`](https://www.cvxgrp.org/CVXR/reference/curvatures.md),
[`is_dgp()`](https://www.cvxgrp.org/CVXR/reference/is_dgp.md)
