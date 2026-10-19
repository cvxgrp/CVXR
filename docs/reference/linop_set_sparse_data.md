# Set sparse data on a LinOp

Set sparse data on a LinOp

## Usage

``` r
linop_set_sparse_data(ptr, value)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- value:

  Sparse matrix data (dgCMatrix or similar).

## Value

No return value (`NULL`), called for its side effect: the sparse data of
the C++ LinOp behind `ptr` is set in place.
