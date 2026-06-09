# Reshape a vector into an upper triangular matrix

Inverts
[`upper_tri`](https://www.cvxgrp.org/CVXR/reference/upper_tri.md). Takes
a flat vector and returns an upper triangular matrix (row-major order,
matching CVXPY convention).

## Usage

``` r
vec_to_upper_tri(expr, strict = FALSE)
```

## Arguments

- expr:

  An Expression (vector).

- strict:

  Logical. If TRUE, returns a strictly upper triangular matrix (diagonal
  is zero). If FALSE, includes the diagonal. Default is FALSE.

## Value

An Expression representing the upper triangular matrix.
