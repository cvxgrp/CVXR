# Check if an Expression is Linearizable-Concave

Concave after linearizing all smooth subexpressions (DNLP composition
rule). Mirrors CVXPY's `Expression.is_linearizable_concave()`.

## Usage

``` r
is_linearizable_concave(x, ...)
```

## Arguments

- x:

  An expression object.

- ...:

  Not used.

## Value

Logical scalar.

## See also

[`is_dnlp`](https://www.cvxgrp.org/CVXR/reference/is_dnlp.md),
[`is_linearizable_convex`](https://www.cvxgrp.org/CVXR/reference/is_linearizable_convex.md)
