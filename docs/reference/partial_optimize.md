# Partial optimization transform

Builds an Expression representing the optimal value of `prob` as a
function of the variables you choose NOT to optimise over. Useful for
two-stage / hierarchical optimisation, custom atom definitions, and
embedding sub-problems inside larger problems.

## Usage

``` r
partial_optimize(
  prob,
  opt_vars = NULL,
  dont_opt_vars = NULL,
  solver = NULL,
  ...
)
```

## Arguments

- prob:

  A [Problem](https://www.cvxgrp.org/CVXR/reference/Problem.md) to
  partially optimise.

- opt_vars:

  Optional list of
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s to
  optimise over.

- dont_opt_vars:

  Optional list of
  [Variable](https://www.cvxgrp.org/CVXR/reference/Variable.md)s to keep
  as free arguments of the resulting expression.

- solver:

  Optional solver name (passed to
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) when the
  PartialProblem is evaluated via
  [`value()`](https://www.cvxgrp.org/CVXR/reference/value.md) or
  [`grad()`](https://www.cvxgrp.org/CVXR/reference/grad.md)).

- ...:

  Additional named arguments forwarded to
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) when
  [`value()`](https://www.cvxgrp.org/CVXR/reference/value.md) /
  [`grad()`](https://www.cvxgrp.org/CVXR/reference/grad.md) are called.

## Value

A `PartialProblem` expression.

## Details

Exactly one of `opt_vars` or `dont_opt_vars` may be `NULL`; the missing
list is taken to be the complement (relative to the full list of
variables in `prob`). If both are supplied, they must together cover
every variable in `prob`.

The returned `PartialProblem` is an Expression with scalar shape: it is
convex when `prob` is DCP with a `Minimize` objective and concave when
DCP with `Maximize`. Embed it like any other expression in a larger
Problem; the larger problem's canonicalizer will pull the inner
objective and constraints into the outer cone form so a single solve
handles both layers.

## See also

[Problem](https://www.cvxgrp.org/CVXR/reference/Problem.md),
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- Variable(3)
t <- Variable(3)
abs_x <- partial_optimize(
  Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t)),
  opt_vars = list(t)
)
## abs_x is now an expression of x alone, equivalent to sum(abs(x)).
} # }
```
