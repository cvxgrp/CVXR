# Get Expression Curvature

Returns the DCP curvature of an expression as a string.

## Usage

``` r
curvature(x)
```

## Arguments

- x:

  A CVXR expression.

## Value

Character: `"CONSTANT"`, `"AFFINE"`, `"CONVEX"`, `"CONCAVE"`, or
`"UNKNOWN"`.

## See also

[`is_convex()`](https://www.cvxgrp.org/CVXR/reference/is_convex.md),
[`is_concave()`](https://www.cvxgrp.org/CVXR/reference/is_concave.md),
[`is_affine()`](https://www.cvxgrp.org/CVXR/reference/is_affine.md),
[`is_constant()`](https://www.cvxgrp.org/CVXR/reference/is_constant.md)
