# Variance of an expression

Computes the variance of an expression, over all entries or along an
axis.

## Usage

``` r
cvxr_var(x, axis = NULL, keepdims = FALSE, ddof = 0)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL (all), 1 (row-wise), or 2 (column-wise).

- keepdims:

  Logical; keep reduced dimension?

- ddof:

  Degrees of freedom correction (default 0, population variance).

## Value

An Expression representing the variance.
