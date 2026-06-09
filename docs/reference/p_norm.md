# General p-norm of an expression

General p-norm of an expression

## Usage

``` r
p_norm(
  x,
  p = 2,
  axis = NULL,
  keepdims = FALSE,
  max_denom = 1024L,
  approx = TRUE
)
```

## Arguments

- x:

  An Expression

- p:

  Numeric exponent (default 2)

- axis:

  NULL (all), 1 (row-wise), or 2 (column-wise)

- keepdims:

  Logical

- max_denom:

  Integer max denominator for rational approx

- approx:

  If TRUE (default), use SOC approximation. If FALSE, use exact power
  cone.

## Value

A Pnorm, PnormApprox, Norm1, or NormInf atom
