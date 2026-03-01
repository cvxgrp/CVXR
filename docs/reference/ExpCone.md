# Create an Exponential Cone Constraint

Constrains \\(x, y, z)\\ to lie in the exponential cone: \$\$K =
\\(x,y,z) \mid y \exp(x/y) \le z,\\ y \> 0\\\$\$

## Usage

``` r
ExpCone(x_expr, y_expr, z_expr, constr_id = NULL)
```

## Arguments

- x_expr:

  A CVXR expression.

- y_expr:

  A CVXR expression.

- z_expr:

  A CVXR expression.

- constr_id:

  Optional integer constraint ID.

## Value

An `ExpCone` constraint object.

## Details

All three arguments must be affine, real, and have the same shape.
