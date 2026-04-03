# Solve a Convex Optimization Problem

Solves the problem and returns the optimal objective value. After
solving, variable values can be retrieved with
[`value`](https://www.cvxgrp.org/CVXR/reference/value.md), constraint
dual values with
[`dual_value`](https://www.cvxgrp.org/CVXR/reference/dual_value.md), and
solver information with
[`solver_stats`](https://www.cvxgrp.org/CVXR/reference/solver_stats.md).

## Usage

``` r
psolve(
  problem,
  solver = NULL,
  gp = FALSE,
  qcp = FALSE,
  verbose = FALSE,
  warm_start = FALSE,
  ...
)
```

## Arguments

- problem:

  A [`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)
  object.

- solver:

  Character string naming the solver to use (e.g., `"CLARABEL"`,
  `"SCS"`, `"OSQP"`, `"HIGHS"`), or `NULL` for automatic selection.

- gp:

  Logical; if `TRUE`, solve as a geometric program (DGP).

- qcp:

  Logical; if `TRUE`, solve as a quasiconvex program (DQCP) via
  bisection. Only needed for non-DCP DQCP problems.

- verbose:

  Logical; if `TRUE`, print solver output.

- warm_start:

  Logical; if `TRUE`, use the current variable values as a warm-start
  point for the solver.

- ...:

  Solver options passed to
  [`solver_opts()`](https://www.cvxgrp.org/CVXR/reference/solver_opts.md).
  Includes chain-construction options (`use_quad_obj`), standard
  tolerances (`feastol`, `reltol`, `abstol`, `num_iter`), and
  solver-specific parameters (e.g., `eps_abs`, `scip_params`). See
  [`solver_opts`](https://www.cvxgrp.org/CVXR/reference/solver_opts.md)
  for details. For DQCP problems (`qcp = TRUE`), additional arguments
  include `low`, `high`, `eps`, `max_iters`, and
  `max_iters_interval_search`.

## Value

The optimal objective value (numeric scalar), or `Inf` / `-Inf` for
infeasible / unbounded problems.

## See also

[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md),
[`status`](https://www.cvxgrp.org/CVXR/reference/status.md),
[`solver_stats`](https://www.cvxgrp.org/CVXR/reference/solver_stats.md),
[`solver_default_param`](https://www.cvxgrp.org/CVXR/reference/solver_default_param.md)

## Examples

``` r
x <- Variable()
prob <- Problem(Minimize(x), list(x >= 5))
result <- psolve(prob, solver = "CLARABEL")
```
