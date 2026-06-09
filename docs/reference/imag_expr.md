# Extract Imaginary Part of Expression

Returns the imaginary part of a complex expression. R's native
[`Im()`](https://rdrr.io/r/base/complex.html) dispatches here for CVXR
expressions via the Complex S3 group handler.

## Usage

``` r
imag_expr(expr)
```

## Arguments

- expr:

  A CVXR Expression.

## Value

An `Imag_` atom.
