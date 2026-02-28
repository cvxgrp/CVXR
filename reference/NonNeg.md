# Create a Non-Negative Constraint

Constrains an expression to be non-negative elementwise: \\x \ge 0\\.

## Usage

``` r
NonNeg(expr, constr_id = NULL)
```

## Arguments

- expr:

  A CVXR expression.

- constr_id:

  Optional integer constraint ID.

## Value

A `NonNeg` constraint object.

## See also

[`NonPos`](https://www.cvxgrp.org/CVXR/reference/NonPos.md),
[`Inequality`](https://www.cvxgrp.org/CVXR/reference/Inequality.md)
