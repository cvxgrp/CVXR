# Get the DCP Sign of an Expression

Returns the sign of an expression under DCP analysis. Use this instead
of [`sign()`](https://rdrr.io/r/base/sign.html), which conflicts with
the base R function.

## Usage

``` r
expr_sign(x, ...)
```

## Arguments

- x:

  An expression object.

- ...:

  Not used.

## Value

Character string: `"POSITIVE"`, `"NEGATIVE"`, `"ZERO"`, or `"UNKNOWN"`.
