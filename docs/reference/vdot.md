# Vector dot product (inner product)

The standard inner product of `x` and `y`: both are flattened,
multiplied elementwise, and summed. **Conjugate-linear in `x`** — `x` is
conjugated before multiplying, matching `numpy.vdot()` and CVXPY's
`vdot()`. For real arguments conjugation is the identity, so this is the
ordinary dot product.

Either argument may also be a (possibly nested) list, which is flattened
before the product; the two are flattened independently, so their
nesting need not match.

## Usage

``` r
vdot(x, y)
```

## Arguments

- x:

  An Expression, numeric value, or nested list thereof. The
  conjugate-linear argument.

- y:

  An Expression, numeric value, or nested list thereof. The linear
  argument.

## Value

A scalar Expression representing `sum(Conj(x) * y)`.

## See also

[`scalar_product()`](https://www.cvxgrp.org/CVXR/reference/scalar_product.md),
[`cvxr_outer()`](https://www.cvxgrp.org/CVXR/reference/cvxr_outer.md),
[`deep_flatten()`](https://www.cvxgrp.org/CVXR/reference/deep_flatten.md)

## Examples

``` r
x <- Variable(3)
vdot(x, c(1, 2, 3))
#> SumEntries(CVXR::Conj_(Reshape(var171, c(3, 1))) * Reshape([3x1 matrix], c(3, 1)), NULL, FALSE)

a <- Variable(); b <- Variable()
vdot(list(a, b), c(1, 2))
#> SumEntries(CVXR::Conj_(CVXR::VStack(Reshape(var178, c(1, 1)), Reshape(var179, c(1, 1)))) * Reshape([2x1 matrix], c(2, 1)), NULL, FALSE)
```
