# Horizontal concatenation of expressions

Named operands are accepted, as in
[`cbind()`](https://rdrr.io/r/base/cbind.html), and the names are
dropped.

## Usage

``` r
hstack(..., deparse.level = 1)
```

## Arguments

- ...:

  Expressions (same number of rows)

- deparse.level:

  Accepted for compatibility with
  [`cbind()`](https://rdrr.io/r/base/cbind.html); expressions have no
  dimnames, so it has no effect

## Value

An HStack atom
