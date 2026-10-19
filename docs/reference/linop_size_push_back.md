# Push a dimension to the LinOp size vector

Push a dimension to the LinOp size vector

## Usage

``` r
linop_size_push_back(ptr, value)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- value:

  Integer dimension value.

## Value

No return value (`NULL`), called for its side effect: `value` is
appended to the size vector of the C++ LinOp behind `ptr`.
