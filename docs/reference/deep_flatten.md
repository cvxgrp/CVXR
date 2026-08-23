# Recursively flatten a nested list of expressions into one column vector

Flattens `x` into a single column vector. An Expression or numeric is
vectorized column-major; a list is flattened element by element and the
pieces stacked in order, recursively. This is what lets
[`vdot()`](https://www.cvxgrp.org/CVXR/reference/vdot.md) accept nested
lists: `vdot(list(a, b), c(1, 2))` is `a * 1 + b * 2`.

## Usage

``` r
deep_flatten(x)
```

## Arguments

- x:

  An Expression, a numeric value, or a (possibly nested) list of them.

## Value

An Expression of shape `c(n, 1)`.

## See also

[`vec()`](https://www.cvxgrp.org/CVXR/reference/vec.md),
[`vdot()`](https://www.cvxgrp.org/CVXR/reference/vdot.md)

## Examples

``` r
a <- Variable(); b <- Variable()
deep_flatten(list(a, b))
#> CVXR::VStack(Reshape(var57, c(1, 1)), Reshape(var58, c(1, 1)))
```
