# Elementwise maximum of expressions

Named operands are accepted, as in
[`pmax()`](https://rdrr.io/r/base/Extremes.html), and the names are
dropped.

## Usage

``` r
max_elemwise(..., na.rm = FALSE)
```

## Arguments

- ...:

  Expressions (at least 2)

- na.rm:

  Accepted for compatibility with
  [`pmax()`](https://rdrr.io/r/base/Extremes.html); expressions cannot
  contain `NA`, so it has no effect

## Value

A Maximum atom
