# Set the data dimensionality on a LinOp

Set the data dimensionality on a LinOp

## Usage

``` r
linop_set_data_ndim(ptr, value)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- value:

  Integer dimensionality (0 = scalar, 2 = matrix).

## Value

No return value (`NULL`), called for its side effect: the data
dimensionality of the C++ LinOp behind `ptr` is set in place.
