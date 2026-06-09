# Elementwise log of the gamma function

Piecewise linear approximation of `log(gamma(x))`. Has modest accuracy
over the full range, approaching perfect accuracy as x goes to infinity.

## Usage

``` r
loggamma(x)
```

## Arguments

- x:

  An Expression or numeric value (must be positive for DCP).

## Value

A convex Expression representing log(gamma(x)).
