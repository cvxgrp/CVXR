# Elementwise log of the standard normal CDF

Quadratic approximation of `log(pnorm(x))` with modest accuracy over the
range -4 to 4.

## Usage

``` r
log_normcdf(x)
```

## Arguments

- x:

  An Expression or numeric value.

## Value

A concave Expression representing log(Phi(x)).
