# Create a Zero Constraint

Constrains an expression to equal zero elementwise: \\x = 0\\.

## Usage

``` r
Zero(expr, constr_id = NULL)
```

## Arguments

- expr:

  A CVXR expression.

- constr_id:

  Optional integer constraint ID.

## Value

A `Zero` constraint object.

## See also

[`Equality`](https://www.cvxgrp.org/CVXR/reference/Equality.md),
[`NonPos`](https://www.cvxgrp.org/CVXR/reference/NonPos.md),
[`NonNeg`](https://www.cvxgrp.org/CVXR/reference/NonNeg.md)
