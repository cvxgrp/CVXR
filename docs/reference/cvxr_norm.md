# Compute a norm of an expression

Compute a norm of an expression

## Usage

``` r
cvxr_norm(x, p = 2, axis = NULL, keepdims = FALSE, max_denom = 1024L)
```

## Arguments

- x:

  An Expression

- p:

  Norm type: 1, 2, Inf, or "fro" (Frobenius)

- axis:

  NULL (all), 0 (columns), or 1 (rows)

- keepdims:

  Logical

- max_denom:

  Integer max denominator

## Value

A norm atom
