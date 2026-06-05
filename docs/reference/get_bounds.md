# Lower/Upper Bounds of a Leaf

Returns the effective `(lower, upper)` bounds of a leaf, combining its
`bounds` attribute with sign (`nonneg`/`pos`/`nonpos`/`neg`) and
`boolean` attributes. Used by the NLP (DNLP) solve path to form variable
bounds.

## Usage

``` r
get_bounds(x, ...)
```

## Arguments

- x:

  An expression (a
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)/leaf, or
  a composite expression whose bounds are propagated from its
  arguments).

- ...:

  Passed to methods; currently unused.

## Value

A list `list(lower, upper)` of two real matrices matching the expression
shape (column-major). For leaves these come from explicit bounds and
sign attributes; for atoms they are propagated from argument bounds via
interval arithmetic (see `bounds_from_args`).
