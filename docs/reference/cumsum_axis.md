# Cumulative sum along an axis

Cumulative sum along an axis

## Usage

``` r
cumsum_axis(x, axis = 2L)
```

## Arguments

- x:

  An Expression

- axis:

  2 (down columns, the default, as CVXPY's `axis = 0`), 1 (across rows),
  or NULL (all entries, in column-major order, as a column)

## Value

A Cumsum atom
