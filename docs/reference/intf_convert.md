# Convert a value to a numeric matrix or sparse matrix

Normalizes R values so the rest of CVXR can assume a consistent type.
Scalars -\> 1x1 matrix, vectors -\> column matrix, logical -\> numeric.
Sparse matrices are kept sparse.

## Usage

``` r
intf_convert(val)
```

## Arguments

- val:

  A numeric scalar, vector, matrix, or Matrix object

## Value

A matrix or dgCMatrix
