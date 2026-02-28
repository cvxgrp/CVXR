# Mean of an expression

Computes the arithmetic mean of an expression along an axis.

## Usage

``` r
cvxr_mean(x, axis = NULL, keepdims = FALSE)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL (all), 1 (row-wise), or 2 (column-wise).

- keepdims:

  Logical; keep reduced dimension?

## Value

An Expression representing the mean.
