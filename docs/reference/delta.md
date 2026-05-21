# Access the perturbation `delta` of a Variable or Parameter

Used by [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md)
with `requires_grad = TRUE` and `Problem$derivative()`. On a
`Parameter`, the user sets `delta` to a perturbation of the parameter's
value;
[`derivative()`](https://www.cvxgrp.org/CVXR/reference/derivative.md)
then reports the predicted change in each `Variable`'s optimal value as
`delta(variable)`.

## Usage

``` r
delta(x)

delta(x) <- value
```

## Arguments

- x:

  A `Variable` or `Parameter`.

- value:

  A numeric array of the same shape as `x`, or `NULL`.

## Value

The perturbation (numeric array) or `NULL`.
