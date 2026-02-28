# Reshape an expression to a new shape

Reshape an expression to a new shape

## Usage

``` r
reshape_expr(x, dim, order = "F")
```

## Arguments

- x:

  An Expression or numeric value.

- dim:

  Integer vector of length 2: the target shape c(nrow, ncol). A single
  integer is treated as c(dim, 1). Use -1 to infer a dimension.

- order:

  Character: "F" (column-major, default) or "C" (row-major).

## Value

A Reshape expression.
