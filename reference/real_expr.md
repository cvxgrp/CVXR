# Extract Real Part of Expression

Returns the real part of a complex expression. R's native
[`Re()`](https://rdrr.io/r/base/complex.html) dispatches here for CVXR
expressions via the Complex S3 group handler.

## Usage

``` r
real_expr(expr)
```

## Arguments

- expr:

  A CVXR Expression.

## Value

A `Real_` atom.
