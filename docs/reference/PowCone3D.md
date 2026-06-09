# Create a 3D Power Cone Constraint

Constrains \\(x, y, z)\\ to lie in the 3D power cone: \$\$x^\alpha \cdot
y^{1-\alpha} \ge \|z\|, \quad x \ge 0, \\ y \ge 0\$\$

## Usage

``` r
PowCone3D(x_expr, y_expr, z_expr, alpha, constr_id = NULL)
```

## Arguments

- x_expr:

  A CVXR expression.

- y_expr:

  A CVXR expression.

- z_expr:

  A CVXR expression.

- alpha:

  A CVXR expression or numeric value in \\(0, 1)\\.

- constr_id:

  Optional integer constraint ID.

## Value

A `PowCone3D` constraint object.

## See also

[`PowConeND`](https://www.cvxgrp.org/CVXR/reference/PowConeND.md)
