# Per-column Subgradient for AxisAtoms (private)

R counterpart of CVXPY's `AxisAtom._column_grad(self, value)`. Receives
a single column (a 1-D fiber) of the input and returns the gradient of
the atom's reduction over that column. Used by `.grad(AxisAtom)`, which
walks the input by fibers and assembles the full Jacobian.

## Usage

``` r
.column_grad(x, value, ...)
```

## Arguments

- x:

  An axis-atom expression.

- value:

  A numeric vector (the fiber).

- ...:

  Reserved for future use; method dispatch ignores it.

## Value

A numeric vector of the same length, or `NULL` if the gradient is
undefined for this fiber.
