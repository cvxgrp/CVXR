# Outer product of two vectors

Computes the outer product `x %*% t(y)`. Both inputs must be vectors.

## Usage

``` r
cvxr_outer(x, y)
```

## Arguments

- x:

  An Expression or numeric value (vector).

- y:

  An Expression or numeric value (vector).

## Value

An Expression of shape (length(x), length(y)).
