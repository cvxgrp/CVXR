# Create an N-Dimensional Power Cone Constraint

Constrains \\(W, z)\\ to lie in the N-dimensional power cone: \$\$\prod
W_i^{\alpha_i} \ge \|z\|, \quad W \ge 0\$\$ where \\\alpha_i \> 0\\ and
\\\sum \alpha_i = 1\\.

## Usage

``` r
PowConeND(W, z, alpha, axis = 2L, constr_id = NULL)
```

## Arguments

- W:

  A CVXR expression (vector or matrix).

- z:

  A CVXR expression (scalar or vector).

- alpha:

  A CVXR expression with positive entries summing to 1 along the
  specified axis.

- axis:

  Integer, 2 (default, column-wise) or 1 (row-wise).

- constr_id:

  Optional integer constraint ID.

## Value

A `PowConeND` constraint object.

## Known limitations

The R clarabel solver does not currently support the PowConeND cone
specification. Problems involving PowConeND (e.g., exact geometric mean
with more than 2 arguments) should use SCS or MOSEK as the solver, or
use approximation-based atoms (e.g., `geo_mean(x, approx = TRUE)`).

## See also

[`PowCone3D`](https://www.cvxgrp.org/CVXR/reference/PowCone3D.md)
