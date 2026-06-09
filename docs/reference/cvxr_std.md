# Standard deviation of an expression

Computes the standard deviation of an expression.

## Usage

``` r
cvxr_std(x, axis = NULL, keepdims = FALSE, ddof = 0)

std(x, axis = NULL, keepdims = FALSE, ddof = 0)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL (all), 1 (row-wise), or 2 (column-wise).

- keepdims:

  Logical; keep reduced dimension?

- ddof:

  Degrees of freedom correction (default 0, population std).

## Value

An Expression representing the standard deviation.
