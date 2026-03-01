# Create an Inequality Constraint

Constrains the left-hand side to be less than or equal to the right-hand
side elementwise: \\lhs \le rhs\\. Typically created via the `<=`
operator on CVXR expressions.

## Usage

``` r
Inequality(lhs, rhs, constr_id = NULL)
```

## Arguments

- lhs:

  A CVXR expression (left-hand side).

- rhs:

  A CVXR expression (right-hand side).

- constr_id:

  Optional integer constraint ID.

## Value

An `Inequality` constraint object.

## See also

[`Equality`](https://www.cvxgrp.org/CVXR/reference/Equality.md),
[`NonPos`](https://www.cvxgrp.org/CVXR/reference/NonPos.md),
[`NonNeg`](https://www.cvxgrp.org/CVXR/reference/NonNeg.md)
