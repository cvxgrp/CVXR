# Access the gradient of a Variable or Parameter

Used by [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md)
with `requires_grad = TRUE` and `Problem$backward()`. Stores a numeric
array of the same shape as the leaf, or `NULL` (the default).

## Usage

``` r
gradient(x)

gradient(x) <- value
```

## Arguments

- x:

  A `Variable` or `Parameter`.

- value:

  A numeric array of the same shape as `x`, or `NULL`.

## Value

The gradient (numeric array) or `NULL`.
