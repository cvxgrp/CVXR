# Get Expression Curvature

Returns the curvature of an expression as a string.

## Usage

``` r
curvature(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

Character: the first curvature that holds, checked in the order
`"CONSTANT"`, `"AFFINE"`, `"CONVEX"`, `"CONCAVE"`, `"LOG-LOG AFFINE"`,
`"LOG-LOG CONVEX"`, `"LOG-LOG CONCAVE"`, `"QUASILINEAR"`,
`"QUASICONVEX"`, `"QUASICONCAVE"`; `"UNKNOWN"` if none does.

## See also

[`curvatures()`](https://www.cvxgrp.org/CVXR/reference/curvatures.md),
[`log_log_curvature()`](https://www.cvxgrp.org/CVXR/reference/log_log_curvature.md),
[`is_convex()`](https://www.cvxgrp.org/CVXR/reference/is_convex.md),
[`is_concave()`](https://www.cvxgrp.org/CVXR/reference/is_concave.md),
[`is_affine()`](https://www.cvxgrp.org/CVXR/reference/is_affine.md),
[`is_constant()`](https://www.cvxgrp.org/CVXR/reference/is_constant.md)
