# Extract Diagonal from a Matrix

Extracts the `k`-th diagonal of a square matrix as a column vector.

## Usage

``` r
DiagMat(x, k = 0L, id = NULL)
```

## Arguments

- x:

  A CVXR expression (square matrix).

- k:

  Integer diagonal offset. `k = 0` (default) is the main diagonal,
  `k > 0` is above, `k < 0` is below.

- id:

  Optional integer ID.

## Value

A `DiagMat` expression of shape `c(n - abs(k), 1)`.

## See also

[`DiagVec`](https://www.cvxgrp.org/CVXR/reference/DiagVec.md)
