# Create a Constant Expression

Wraps a numeric value as a CVXR constant for use in optimization
expressions. Constants are typically created implicitly when combining
numeric values with CVXR expressions via arithmetic operators.

## Usage

``` r
Constant(value, name = NULL)
```

## Arguments

- value:

  A numeric scalar, vector, matrix, or sparse matrix.

- name:

  Optional character string name.

## Value

A `Constant` object (inherits from `Leaf` and `Expression`).

## Examples

``` r
c1 <- Constant(5)
c2 <- Constant(matrix(1:6, 2, 3))
```
