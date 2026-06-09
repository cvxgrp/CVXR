# Positive Semidefinite Constraint Operator

Creates a PSD constraint: `e1 - e2` is positive semidefinite. This is
the R equivalent of Python's `A >> B`.

## Usage

``` r
e1 %>>% e2
```

## Arguments

- e1, e2:

  CVXR expressions or numeric matrices.

## Value

A `PSD` constraint object.

## See also

[`PSD()`](https://www.cvxgrp.org/CVXR/reference/PSD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- Variable(3, 3, symmetric = TRUE)
constr <- X %>>% diag(3)  # X - I is PSD
} # }
```
