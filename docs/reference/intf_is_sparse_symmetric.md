# Check if a sparse matrix is symmetric (or Hermitian)

Check if a sparse matrix is symmetric (or Hermitian)

## Usage

``` r
intf_is_sparse_symmetric(m, complex = FALSE)
```

## Arguments

- m:

  A square sparse matrix

- complex:

  If `TRUE`, test `m == Conj(t(m))`

## Value

Logical scalar
