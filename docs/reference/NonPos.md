# Create a Non-Positive Constraint

Constrains an expression to be non-positive elementwise: \\x \le 0\\.

## Usage

``` r
NonPos(expr, constr_id = NULL)
```

## Arguments

- expr:

  A CVXR expression.

- constr_id:

  Optional integer constraint ID.

## Value

A `NonPos` constraint object.

## See also

[`NonNeg`](https://www.cvxgrp.org/CVXR/reference/NonNeg.md),
[`Inequality`](https://www.cvxgrp.org/CVXR/reference/Inequality.md)
