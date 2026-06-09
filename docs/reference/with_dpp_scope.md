# Execute Expression Within DPP Scope

Within this scope, `Parameter` objects are treated as affine (not
constant) for curvature analysis. This is used internally by
[`is_dpp()`](https://www.cvxgrp.org/CVXR/reference/is_dpp.md) to check
DPP compliance.

## Usage

``` r
with_dpp_scope(expr)
```

## Arguments

- expr:

  An R expression to evaluate within the DPP scope.

## Value

The result of evaluating `expr`.
