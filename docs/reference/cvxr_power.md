# Create a Power atom

CVXPY's `power(x, p)`: `x` raised to the power `p`, for any argument (a
plain number becomes a constant). `cvxr_power()` is the preferred name.
`power()` is the same function, kept because existing code uses it; it
masks [`power`](https://rdrr.io/r/stats/power.html), the unrelated
[`glm()`](https://rdrr.io/r/stats/glm.html) link constructor, so with
CVXR attached write
[`stats::power()`](https://rdrr.io/r/stats/power.html) in a
[`glm()`](https://rdrr.io/r/stats/glm.html) family. `x^p` is the same
atom for an expression `x`.

## Usage

``` r
cvxr_power(x, p, max_denom = 1024L, approx = TRUE)

power(x, p, max_denom = 1024L, approx = TRUE)
```

## Arguments

- x:

  An Expression (the base), OR a positive constant if `p` is a variable
  (the identity `b^x = exp(x * log(b))` is used).

- p:

  Numeric exponent, Parameter, or Expression. If `p` is a non-constant
  Expression and `x` is a positive constant, dispatches to
  `exp(p * log(x))`.

- max_denom:

  Maximum denominator for rational approximation

- approx:

  If TRUE (default), use SOC approximation. If FALSE, use exact power
  cone.

## Value

A Power or PowerApprox atom, or an exp expression for the const-base
case.

## Note

`sqrt(x)` on a CVXR expression dispatches to `Power(x, 0.5)` via the
Math group generic. See
[`math_atoms`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) for
all standard R function dispatch.
