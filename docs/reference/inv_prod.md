# Reciprocal of product of entries

Computes the reciprocal of the product of entries. Equivalent to
`geo_mean(x)^(-n)` where n is the number of entries.

## Usage

``` r
inv_prod(x, approx = TRUE)
```

## Arguments

- x:

  An Expression or numeric value (must have positive entries).

- approx:

  Logical; if TRUE (default), use SOC approximation.

## Value

A convex Expression representing the reciprocal product.
