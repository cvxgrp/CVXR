# Set dense data on a LinOp

Set dense data on a LinOp

## Usage

``` r
linop_set_dense_data(ptr, value)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- value:

  Dense matrix data.

## Value

No return value (`NULL`), called for its side effect: the dense data of
the C++ LinOp behind `ptr` is set in place.
