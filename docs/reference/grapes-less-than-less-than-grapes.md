# Negative Semidefinite Constraint Operator

Creates an NSD constraint: `e2 - e1` is positive semidefinite, i.e.,
`e1` is NSD relative to `e2`. This is the R equivalent of Python's
`A << B`.

## Usage

``` r
e1 %<<% e2
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
constr <- X %<<% diag(3)  # I - X is PSD (X is NSD relative to I)
} # }
```
