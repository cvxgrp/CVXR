# Variance of an expression

Computes the variance. Only supports full reduction (axis = NULL).

## Usage

``` r
cvxr_var(x, axis = NULL, keepdims = FALSE, ddof = 0)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL only (axis reduction not yet supported).

- keepdims:

  Logical; keep reduced dimension?

- ddof:

  Degrees of freedom correction (default 0, population variance).

## Value

An Expression representing the variance.
