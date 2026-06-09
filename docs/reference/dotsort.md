# Weighted sorted dot product

Computes `<sort(vec(X)), sort(vec(W))>` where X is an expression and W
is a constant. A generalization of `sum_largest` and `sum_smallest`.

## Usage

``` r
dotsort(X, W)
```

## Arguments

- X:

  An Expression or numeric value.

- W:

  A constant numeric vector or matrix.

## Value

A scalar convex Expression.
