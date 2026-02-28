# Elementwise Complex Conjugate

Returns the complex conjugate of an expression. For real expressions
this is a no-op. R's native
[`Conj()`](https://rdrr.io/r/base/complex.html) dispatches here for CVXR
expressions via the Complex S3 group handler.

## Usage

``` r
conj_expr(expr)
```

## Arguments

- expr:

  A CVXR Expression.

## Value

A `Conj_` atom.
