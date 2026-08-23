# Outer product of two vectors

Computes the outer product `x %*% t(y)`. Both inputs must be vectors.
Named `cvxr_outer()` because base R already exports
[`outer()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md).

## Usage

``` r
cvxr_outer(x, y)
```

## Arguments

- x:

  An Expression or numeric value (vector).

- y:

  An Expression or numeric value (vector).

## Value

An Expression of shape (length(x), length(y)).

## See also

[`vdot()`](https://www.cvxgrp.org/CVXR/reference/vdot.md)
