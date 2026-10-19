# Convert a value to a numeric matrix or sparse matrix

Normalizes R values so the rest of CVXR can assume a consistent type, as
CVXPY's
[`Constant()`](https://www.cvxgrp.org/CVXR/reference/Constant.md) does
through `const_to_matrix` (constant.py:51-54;
numpy_interface/sparse_matrix_interface.py:40-51, ndarray_interface.py):
every sparse input becomes a double CSC matrix, bool and int are cast to
double, the shape is kept and explicit zeros are kept.

## Usage

``` r
intf_convert(val)
```

## Arguments

- val:

  A logical, integer, double or complex scalar, vector or matrix, or a
  Matrix object

## Value

A base matrix or a dgCMatrix

## Details

Here any `sparseMatrix` or `sparseVector` becomes a general double
`dgCMatrix` (symmetric, triangular, diagonal, triplet, row-compressed,
pattern and logical classes included); a dense `Matrix` becomes a base
matrix; a base logical, integer, double or complex value is stored as
double (complex kept), a scalar or vector as a column.
