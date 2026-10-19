# Compute a norm of an expression

With `axis = NULL`, `p = 1`, `2` and `Inf` are the matrix norms, as in
CVXPY and base R
[`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md): the
maximum absolute column sum, the largest singular value
([`sigma_max`](https://www.cvxgrp.org/CVXR/reference/sigma_max.md)) and
the maximum absolute row sum. A row vector is a 1 x n matrix, so its
1-norm is its largest absolute entry and its Inf-norm the sum of
absolute entries. For a column vector these equal the entrywise norms.
Any other `p` is an error for a matrix or row vector, as in CVXPY
("Unsupported matrix norm."); use
[`p_norm`](https://www.cvxgrp.org/CVXR/reference/p_norm.md) for the
entrywise p-norm. On a column vector, and with an `axis`, it is the
entrywise norm. The masking
[`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) gives
the same results under base R's type names.

## Usage

``` r
cvxr_norm(x, p = 2, axis = NULL, keepdims = FALSE, max_denom = 1024L)
```

## Arguments

- x:

  An Expression

- p:

  Norm type: 1, 2, Inf, "fro" (Frobenius), or "nuc" (nuclear norm of a
  matrix). "fro" and "nuc" are whole-matrix norms and take no `axis`.

- axis:

  NULL (all), 0 (columns), or 1 (rows)

- keepdims:

  Logical

- max_denom:

  Integer max denominator

## Value

A norm atom
