# Create an Equality Constraint

Constrains two expressions to be equal elementwise: \\lhs = rhs\\.
Typically created via the `==` operator on CVXR expressions.

## Usage

``` r
Equality(lhs, rhs, constr_id = NULL)
```

## Arguments

- lhs:

  A CVXR expression (left-hand side).

- rhs:

  A CVXR expression (right-hand side).

- constr_id:

  Optional integer constraint ID.

## Value

An `Equality` constraint object.

## See also

[`Zero`](https://www.cvxgrp.org/CVXR/reference/Zero.md),
[`Inequality`](https://www.cvxgrp.org/CVXR/reference/Inequality.md)
