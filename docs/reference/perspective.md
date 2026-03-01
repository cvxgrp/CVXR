# Perspective Transform

Creates the perspective transform of a scalar convex or concave
expression. Given a scalar expression `f(x)` and a nonneg variable `s`,
the perspective is `s * f(x/s)`.

## Usage

``` r
perspective(f, s, f_recession = NULL)
```

## Arguments

- f:

  A scalar convex or concave Expression.

- s:

  A nonneg Variable (scalar).

- f_recession:

  Optional recession function for handling `s = 0`.

## Value

A `Perspective` expression.
