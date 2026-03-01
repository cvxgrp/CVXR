# Convert a value to a CVXR Expression

Wraps numeric vectors, matrices, and Matrix package objects as CVXR
[Constant](https://www.cvxgrp.org/CVXR/reference/Constant.md) objects.
Values that are already CVXR expressions are returned unchanged. This is
useful when you need to ensure that one operand in an arithmetic
expression is a CVXR type so that CVXR's operator dispatch takes effect
— for example, when multiplying a sparse Matrix by a Variable with
`%*%`.

## Usage

``` r
as_cvxr_expr(x)
```

## Arguments

- x:

  A numeric vector, matrix,
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix-class.html)
  object, or CVXR expression.

## Value

A CVXR expression (either the input unchanged or wrapped in
[Constant](https://www.cvxgrp.org/CVXR/reference/Constant.md)).

## Examples

``` r
x <- Variable(3)
## Sparse Matrix %*% Variable needs one CVXR operand for dispatch:
A <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1.0)
expr <- as_cvxr_expr(A) %*% x
```
