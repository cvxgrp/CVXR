# Push a slice (index vector) to a LinOp

Push a slice (index vector) to a LinOp

## Usage

``` r
linop_slice_push_back(ptr, int_vector)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- int_vector:

  Integer vector of indices (0-based for C++).

## Value

No return value (`NULL`), called for its side effect: `int_vector` is
appended to the slice list of the C++ LinOp behind `ptr`.
