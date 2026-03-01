# FiniteSet Constraint

Constrain each entry of an Expression to take a value in a given finite
set of real numbers.

## Usage

``` r
FiniteSet(expre, vec, ineq_form = FALSE, constr_id = NULL)
```

## Arguments

- expre:

  An affine Expression.

- vec:

  A numeric vector (or set) of allowed values.

- ineq_form:

  Logical; controls MIP canonicalization strategy. If `FALSE` (default),
  uses equality formulation (one-hot binary). If `TRUE`, uses inequality
  formulation (sorted differences + ordering).

- constr_id:

  Optional integer constraint ID (internal use).

## Value

A `FiniteSet` constraint.
