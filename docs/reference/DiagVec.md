# Vector to Diagonal Matrix

Constructs a diagonal matrix from a column vector. If `k != 0`, the
vector is placed on the `k`-th super- or sub-diagonal.

## Usage

``` r
DiagVec(x, k = 0L, id = NULL)
```

## Arguments

- x:

  A CVXR expression (column vector).

- k:

  Integer diagonal offset. `k = 0` (default) is the main diagonal,
  `k > 0` is above, `k < 0` is below.

- id:

  Optional integer ID.

## Value

A `DiagVec` expression of shape `c(n + abs(k), n + abs(k))`.

## See also

[`DiagMat`](https://www.cvxgrp.org/CVXR/reference/DiagMat.md)
