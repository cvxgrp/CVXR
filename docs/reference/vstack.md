# Vertical concatenation of expressions

Named operands are accepted, as in
[`rbind()`](https://rdrr.io/r/base/cbind.html), and the names are
dropped.

## Usage

``` r
vstack(..., deparse.level = 1)
```

## Arguments

- ...:

  Expressions (same number of columns)

- deparse.level:

  Accepted for compatibility with
  [`rbind()`](https://rdrr.io/r/base/cbind.html); expressions have no
  dimnames, so it has no effect

## Value

A VStack atom
