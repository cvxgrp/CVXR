# Convert a value to a CVXR Expression

Wraps numeric vectors, matrices, and Matrix package objects as CVXR
[Constant](https://www.cvxgrp.org/CVXR/reference/Constant.md) objects.
Values that are already CVXR expressions are returned unchanged.

## Usage

``` r
as_cvxr_expr(x)
```

## Arguments

- x:

  A numeric vector, matrix,
  [Matrix::Matrix](https://rdrr.io/pkg/Matrix/man/Matrix-class.html)
  object,
  [Matrix::sparseVector](https://rdrr.io/pkg/Matrix/man/sparseVector-class.html)
  object, or CVXR expression.

## Value

A CVXR expression (either the input unchanged or wrapped in
[Constant](https://www.cvxgrp.org/CVXR/reference/Constant.md)).

## Matrix package interoperability

Objects from the Matrix package (`dgCMatrix`, `dgeMatrix`, `ddiMatrix`,
`sparseVector`, etc.) are S4 classes. Because S4 dispatch preempts S7/S3
dispatch, **raw Matrix objects cannot be used directly with CVXR
operators** (`+`, `-`, `*`, `/`, `%*%`, `>=`, `==`, etc.).

Use `as_cvxr_expr()` to wrap a Matrix object as a CVXR
[Constant](https://www.cvxgrp.org/CVXR/reference/Constant.md) before
combining it with CVXR variables or expressions. This preserves sparsity
(unlike [`as.matrix()`](https://rdrr.io/r/base/matrix.html), which
densifies).

Base R `matrix` and `numeric` objects work natively with CVXR operators
— no wrapping is needed.

## Examples

``` r
x <- Variable(3)

## Sparse Matrix needs as_cvxr_expr() for CVXR operator dispatch:
A <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1.0)
expr <- as_cvxr_expr(A) %*% x

## All operators work with wrapped Matrix objects:
y <- Variable(c(3, 3))
expr2 <- as_cvxr_expr(A) + y
constr <- as_cvxr_expr(A) >= y

## Base R matrix works natively (no wrapping needed):
D <- matrix(1:9, 3, 3)
expr3 <- D %*% x
```
