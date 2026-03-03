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
  feastol = NULL,
  reltol = NULL,
  abstol = NULL,
  num_iter = NULL,
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

- feastol:

  Numeric or `NULL`; feasibility tolerance. Mapped to the
  solver-specific parameter name (e.g., `tol_feas` for Clarabel,
  `eps_prim_inf` for OSQP). If `NULL` (default), the solver's own
  default is used.

- reltol:

  Numeric or `NULL`; relative tolerance. Mapped to the solver-specific
  parameter name (e.g., `tol_gap_rel` for Clarabel, `eps_rel` for
  OSQP/SCS). If `NULL` (default), the solver's own default is used.

- abstol:

  Numeric or `NULL`; absolute tolerance. Mapped to the solver-specific
  parameter name (e.g., `tol_gap_abs` for Clarabel, `eps_abs` for
  OSQP/SCS). If `NULL` (default), the solver's own default is used.

- num_iter:

  Integer or `NULL`; maximum number of solver iterations. Mapped to the
  solver-specific parameter name (e.g., `max_iter` for Clarabel/OSQP,
  `max_iters` for SCS). If `NULL` (default), the solver's own default is
  used.

- ...:

  Additional solver-specific options passed directly to the solver. If a
  solver-native parameter name conflicts with a standard parameter
  (e.g., both `feastol = 1e-3` and `tol_feas = 1e-6` are given), the
  solver-native name in `...` takes priority. For DQCP problems
  (`qcp = TRUE`), additional arguments include `low`, `high`, `eps`,
  `max_iters`, and `max_iters_interval_search`.

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
