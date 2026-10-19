# Set a LinOp data sub-tree on a LinOp

Set a LinOp data sub-tree on a LinOp

## Usage

``` r
linop_set_linop_data(ptr, data_ptr)
```

## Arguments

- ptr:

  External pointer to the LinOp.

- data_ptr:

  External pointer to the data LinOp sub-tree.

## Value

No return value (`NULL`), called for its side effect: the sub-tree
becomes the data of the C++ LinOp behind `ptr`.
