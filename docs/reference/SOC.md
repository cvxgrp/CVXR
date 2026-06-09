# Create a Second-Order Cone Constraint

Constrains \\\\X_i\\\_2 \le t_i\\ for each column or row \\i\\, where
`t` is a vector and `X` is a matrix.

## Usage

``` r
SOC(t, X, axis = 2L, constr_id = NULL)
```

## Arguments

- t:

  A CVXR expression (scalar or vector) representing the upper bound.

- X:

  A CVXR expression (vector or matrix) whose columns/rows are bounded.

- axis:

  Integer, 2 (default, column-wise) or 1 (row-wise). Determines whether
  columns (`2`) or rows (`1`) of `X` define the individual cones.

- constr_id:

  Optional integer constraint ID.

## Value

An `SOC` constraint object.
