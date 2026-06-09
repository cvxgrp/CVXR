# Conjugate-Transpose of an Expression

Equivalent to CVXPY's `.H` property. For real expressions, returns
`t(x)`. For complex expressions, returns `Conj(t(x))`.

## Usage

``` r
expr_H(x)
```

## Arguments

- x:

  A CVXR Expression.

## Value

The conjugate-transpose expression.
