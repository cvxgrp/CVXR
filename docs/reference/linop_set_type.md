# Set the type of a LinOp

Set the type of a LinOp

## Usage

``` r
linop_set_type(ptr, type)
```

## Arguments

- ptr:

  External pointer to a LinOp.

- type:

  Character string – one of `.LINOP_TYPES`.

## Value

No return value (`NULL`), called for its side effect: the type of the
C++ LinOp behind `ptr` is set in place.
