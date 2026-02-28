# Create a Power atom

Create a Power atom

## Usage

``` r
power(x, p, max_denom = 1024L, approx = TRUE)
```

## Arguments

- x:

  An Expression

- p:

  Numeric exponent

- max_denom:

  Maximum denominator for rational approximation

- approx:

  If TRUE (default), use SOC approximation. If FALSE, use exact power
  cone.

## Value

A Power or PowerApprox atom

## Note

`sqrt(x)` on a CVXR expression dispatches to `Power(x, 0.5)` via the
Math group generic. See
[`math_atoms`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) for
all standard R function dispatch.
