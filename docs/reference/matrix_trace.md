# Trace of a square matrix expression

For `matrix_trace(A %*% B)`, uses the O(n^2) identity
`trace(A %*% B) = sum(A * t(B))` instead of forming the full matrix
product.

## Usage

``` r
matrix_trace(x)
```

## Arguments

- x:

  An Expression (square matrix)

## Value

A Trace atom or equivalent expression (scalar)
