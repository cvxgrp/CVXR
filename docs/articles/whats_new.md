# What's New in CVXR

This vignette highlights notable user-facing changes in each release of
the S7 rewrite of `CVXR`, newest first. For the complete, fine-grained
list see the package `NEWS` file (`news(package = "CVXR")`). For the
introductory tutorial, see
[`vignette("cvxr_intro")`](https://www.cvxgrp.org/CVXR/articles/cvxr_intro.md);
for worked examples, visit the [CVXR website](https://cvxr.rbind.io).

- [CVXR 1.9.1](#cvxr-19x) — disciplined nonlinear programming,
  derivatives, bounds propagation, new atoms
- [CVXR 1.8.x](#cvxr-18x) — the ground-up S7 rewrite

## CVXR 1.9.1

`CVXR` 1.9.1 is the first CRAN release since 1.8.2 and is a large one:
it folds in the internal 1.8.2-1 and 1.9.0 development cycles. Its
headline additions are disciplined nonlinear programming, a derivative /
sensitivity-analysis API, and interval-bounds propagation with native
solver-bound support.

### Disciplined Nonlinear Programming (DNLP)

`CVXR` 1.9.1 extends modeling beyond convex optimization to **smooth
nonlinear programs**, which need not be convex. You build the problem
from differentiable atoms, check it with
[`is_dnlp()`](https://www.cvxgrp.org/CVXR/reference/is_dnlp.md), and
solve it with `psolve(prob, nlp = TRUE)`. Every DCP problem is also a
DNLP, and the disciplined nonlinear grammar additionally allows smooth
atoms in forms DCP forbids (for example, a product of two
variable-dependent expressions).

`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``2``)`` ``prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(`[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum_squares`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)`(``x`` ``-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2``)``)``)``)`` `[`is_dnlp`](https://www.cvxgrp.org/CVXR/reference/is_dnlp.md)`(``prob``)`` ``# TRUE`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``, nlp ``=`` ``TRUE``)`` ``# solved through the NLP path`

- **New smooth atoms** usable anywhere in a DNLP:
  [`sin()`](https://rdrr.io/r/base/Trig.html),
  [`cos()`](https://rdrr.io/r/base/Trig.html),
  [`tan()`](https://rdrr.io/r/base/Trig.html),
  [`sinh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`asinh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`atanh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`normcdf()`](https://www.cvxgrp.org/CVXR/reference/normcdf.md), and
  [`prod()`](https://rdrr.io/r/base/prod.html).
- **Nonconvex problems** may have several local optima. `best_of = n`
  solves from `n` random initial points (drawn from variable bounds or
  [`sample_bounds()`](https://www.cvxgrp.org/CVXR/reference/sample_bounds.md))
  and keeps the best.
- The NLP path is powered by automatic derivatives from the optional
  `sparsediff` package and a nonlinear solver. Two are supported, both
  in `Enhances` (so guard their use with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)): **UNO**
  (the `Uno` package — headed for CRAN, self-contained, and the
  practical default for R users) and **IPOPT** (the `ipopt` package —
  not on CRAN owing to licensing, so rarely installed). With
  `nlp = TRUE`, IPOPT is preferred when present (matching CVXPY),
  otherwise UNO is used. The `UNO` path also recovers constraint duals
  via
  [`dual_value()`](https://www.cvxgrp.org/CVXR/reference/dual_value.md);
  the `IPOPT` path returns none, matching CVXPY.

See the [DNLP
Tutorial](https://cvxr.rbind.io/examples/dnlp/tutorial.html) for worked
examples.

### Derivatives and sensitivity analysis

`CVXR` 1.9.1 adds the ability to **differentiate the solution map** of a
disciplined problem — to see how the optimal solution responds to small
changes in the parameters (sensitivity analysis) and to compute
gradients of scalar functions of the solution. Request derivatives at
solve time with `requires_grad = TRUE`.

**Forward mode** (perturb parameters, see the change in the solution):

[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, requires_grad ``=`` ``TRUE``)`` `[`delta`](https://www.cvxgrp.org/CVXR/reference/delta.md)`(``a``)`` ``<-`` ``da`` ``# perturbation of parameter a`` `[`derivative`](https://www.cvxgrp.org/CVXR/reference/derivative.md)`(``problem``)`` ``# propagate forward`` `[`delta`](https://www.cvxgrp.org/CVXR/reference/delta.md)`(``x``)`` ``# resulting change in variable x`

**Reverse mode** (gradient of the solution with respect to parameters):

[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, requires_grad ``=`` ``TRUE``)`` `[`backward`](https://www.cvxgrp.org/CVXR/reference/backward.md)`(``problem``)`` ``# propagate backward`` `[`gradient`](https://www.cvxgrp.org/CVXR/reference/gradient.md)`(``a``)`` ``# d(solution) / d(a)`

The chain rule is wired through the `Dgp2Dcp` (log/exp) and
`Complex2Real` reductions, so geometric and complex problems
differentiate too. The derivative API is backed by the optional `diffcp`
R package. See the
[Derivatives](https://cvxr.rbind.io/examples/derivatives/fundamentals.html)
examples and [Sensitivity
Analysis](https://cvxr.rbind.io/examples/dpp/sensitivity-analysis.html).

### Bounds propagation and richer variable bounds

[`get_bounds()`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)
now works on **any expression**, not just variables, propagating
interval bounds through affine, elementwise, and piecewise-linear atoms:

`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``3``, bounds ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``-``1``, ``2``)``)`` `[`get_bounds`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)`(``A`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``x`` ``+`` ``b``)`` ``# bounds propagated through the affine map`` `[`get_bounds`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(``x``)``)`` ``# and through atoms`

Variable bounds may also be **sparse `Matrix` objects** or **symbolic
bounds** involving `Parameter`s; symbolic bounds are enforced at solve
time and update on DPP re-solves. Positive (DGP) variables accept
numeric *and* parametric bounds under `gp = TRUE`.

### New atoms and DPP refinements

- **[`convolve()`](https://www.cvxgrp.org/CVXR/reference/convolve.md)**
  — the (numpy-style) name for the 1-D discrete-convolution
  [`conv()`](https://www.cvxgrp.org/CVXR/reference/conv.md) atom; falls
  through to
  [`stats::convolve()`](https://rdrr.io/r/stats/convolve.html) on
  numeric input.
- **[`is_dpp()`](https://www.cvxgrp.org/CVXR/reference/is_dpp.md) gains
  a `context` argument** (`"dcp"` or `"dgp"`), matching CVXPY’s
  `is_dpp(context = ...)`.

### Solvers

- **CPLEX** now solves LP / SOCP / MI-LP / MI-SOCP through the conic
  path.
- **Native variable bounds**: HiGHS, Gurobi, CPLEX, XPRESS, PIQP, and
  SCIP now consume dense numeric variable bounds directly (including
  parametric bounds for HiGHS), avoiding extra bound constraints and
  speeding up DPP re-solves.

### Bug fixes (also affected 1.8.x)

- [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
  /
  [`get_problem_data()`](https://www.cvxgrp.org/CVXR/reference/get_problem_data.md)
  now take an explicit `gp` argument; previously `gp = TRUE` passed
  through `...` was silently ignored, compiling a geometric program as a
  DCP problem.
- [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) now
  takes explicit `enforce_dpp` and `ignore_dpp` arguments, matching
  CVXPY’s [`solve()`](https://rdrr.io/r/base/solve.html); previously
  they were silently swallowed by `...`.

### Performance

Canonicalization and solving are now faster than the 1.8.2 CRAN release
(roughly 5–13% lower wall-clock on solve-dominated problems such as many
small constraints, SOCPs, and Kalman smoothing), with deterministic
memory allocation unchanged.

## CVXR 1.8.x

### Complete Rewrite Using S7

CVXR 1.8.x is a ground-up rewrite using R’s
[S7](https://rconsortium.github.io/S7/) object system, designed to be
isomorphic with [CVXPY 1.8.2](https://www.cvxpy.org/) for long-term
maintainability. It is approximately 4–5x faster than the previous
S4-based release. This section summarizes the key changes from CVXR 1.x
that may affect users.

### New Features

- **S7 class system** replaces S4 for all expression, constraint, and
  problem classes. Significantly faster construction and method
  dispatch.
- **15 solvers**: CLARABEL (default), SCS, OSQP, HiGHS, MOSEK, Gurobi,
  GLPK, GLPK_MI, ECOS, ECOS_BB, CPLEX, CVXOPT, PIQP, SCIP, and XPRESS.
- **Mixed-integer programming** via GLPK_MI, ECOS_BB, Gurobi, CPLEX,
  HiGHS, SCIP, or XPRESS (`boolean = TRUE` or `integer = TRUE` in
  [`Variable()`](https://www.cvxgrp.org/CVXR/reference/Variable.md)).
- **Parameter support** via
  [`Parameter()`](https://www.cvxgrp.org/CVXR/reference/Parameter.md)
  class and `EvalParams` reduction.
- **50+ atom classes** covering LP, QP, SOCP, SDP, exponential cone, and
  power cone problems.
- **DPP** (Disciplined Parameterized Programming) for efficient
  parameter re-solve with compilation caching.
- **DGP** (Disciplined Geometric Programming) via
  `psolve(prob, gp = TRUE)`.
- **DQCP** (Disciplined Quasiconvex Programming) via
  `psolve(prob, qcp = TRUE)`.
- **Complex variable support** via `Variable(n, complex = TRUE)`.
- **Warm-start support** for several solvers (OSQP, SCS, Gurobi, MOSEK,
  CLARABEL, HiGHS).
- **Matrix package interoperability** via
  [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md).
  Matrix package objects (`dgCMatrix`, `dgeMatrix`, `dsCMatrix`,
  `ddiMatrix`, `sparseVector`) use S4 dispatch which preempts S7/S3, so
  they cannot be used directly with CVXR operators. Wrapping with
  [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md)
  converts them to CVXR `Constant` objects while preserving sparsity
  (unlike [`as.matrix()`](https://rdrr.io/r/base/matrix.html) which
  densifies). Base R `matrix` and `numeric` objects work natively
  without wrapping.

### New solve interface

The primary solve function is now
[`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md), which
returns the optimal value directly:

[`library`](https://rdrr.io/r/base/library.html)`(`[`CVXR`](https://cvxr.rbind.io)`)`` ``x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``2``)`` ``prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(`[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum_squares`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)`(``x``)``)``, `[`list`](https://rdrr.io/r/base/list.html)`(``x`` ``>=`` ``1``)``)`` ``opt_val`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``)`` ``# returns optimal value directly`` ``x_val`` ``<-`` `[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``x``)`` ``# extract variable value`` ``prob_status`` ``<-`` `[`status`](https://www.cvxgrp.org/CVXR/reference/status.md)`(``prob``)`` ``# check status`

The old [`solve()`](https://rdrr.io/r/base/solve.html) still works but
returns a backward-compatible list:

`result`` ``<-`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``prob``)`` ``result``$``value`` ``# optimal value`` ``result``$``getValue``(``x``)`` ``# variable value (deprecated)`` ``result``$``status`` ``# problem status`

### Breaking Changes from CVXR 1.x

#### API changes

| Old API | New API |
|----|----|
| `solve(problem)` | `psolve(problem)` |
| `result$getValue(x)` | `value(x)` |
| `result$value` | return value of [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) |
| `result$status` | `status(problem)` |
| `result$getDualValue(con)` | `dual_value(con)` |
| `problem_status(prob)` | `status(prob)` |
| `problem_solution(prob)` | `solution(prob)` |
| `get_problem_data(prob, solver)` | `problem_data(prob, solver)` |

#### Axis parameter changes

The `axis` parameter now uses R’s
[`apply()`](https://rdrr.io/r/base/apply.html) convention (1-based
indexing):

| Old CVXR    | New CVXR      | Meaning                           |
|-------------|---------------|-----------------------------------|
| `axis = 1`  | `axis = 1`    | Row-wise reduction (unchanged)    |
| `axis = 2`  | `axis = 2`    | Column-wise reduction (unchanged) |
| `axis = NA` | `axis = NULL` | All entries                       |

Passing `axis = 0` now produces an informative error with migration
guidance.

#### PSD constraints

PSD constraints use `PSD(A - B)` instead of `A %>>% B` (though `%>>%`
and `%<<%` operators are still available for backward compatibility).

#### Solver changes

- **Removed**: CBC
- **Added**: HiGHS (LP, QP, MILP), Gurobi (LP, QP, SOCP, MIP), CVXOPT
  (LP, SOCP), PIQP (QP), SCIP, and XPRESS
- **Default solver**: CLARABEL (replaces ECOS)

#### Supported solvers

| Solver   | R Package   | Type     | Problem Classes                     |
|----------|-------------|----------|-------------------------------------|
| CLARABEL | `clarabel`  | Conic    | LP, QP, SOCP, SDP, ExpCone, PowCone |
| SCS      | `scs`       | Conic    | LP, QP, SOCP, SDP, ExpCone, PowCone |
| MOSEK    | `Rmosek`    | Conic    | LP, QP, SOCP, SDP, ExpCone, PowCone |
| ECOS     | `ECOSolveR` | Conic    | LP, SOCP, ExpCone                   |
| ECOS_BB  | `ECOSolveR` | Conic    | LP, SOCP, ExpCone + MI              |
| GUROBI   | `gurobi`    | Conic/QP | LP, QP, SOCP, MI                    |
| GLPK     | `Rglpk`     | Conic    | LP                                  |
| GLPK_MI  | `Rglpk`     | Conic    | LP, MILP                            |
| HIGHS    | `highs`     | Conic/QP | LP, QP, MILP                        |
| CVXOPT   | `cccp`      | Conic    | LP, SOCP                            |
| OSQP     | `osqp`      | QP       | LP, QP                              |
| CPLEX    | `Rcplex`    | Conic/QP | LP, QP, SOCP, MI                    |
| PIQP     | `piqp`      | QP       | LP, QP                              |
| SCIP     | `scip`      | Conic    | LP, MILP, SOCP, MI-SOCP             |
| XPRESS   | `xpress`    | Conic/QP | LP, QP, SOCP, MI                    |

Smooth nonlinear programs additionally use the `IPOPT` and `UNO` NLP
solvers (see [CVXR 1.9.1](#cvxr-19x)).

### New Atoms and Functions

#### Convenience atoms

| Function           | Description                             |
|--------------------|-----------------------------------------|
| `ptp(x)`           | Peak-to-peak (range): `max(x) - min(x)` |
| `cvxr_mean(x)`     | Arithmetic mean along an axis           |
| `cvxr_std(x)`      | Standard deviation                      |
| `cvxr_var(x)`      | Variance                                |
| `vdot(x, y)`       | Vector dot product (inner product)      |
| `cvxr_outer(x, y)` | Outer product of two vectors            |
| `inv_prod(x)`      | Reciprocal of product of entries        |
| `loggamma(x)`      | Elementwise log of gamma function       |
| `log_normcdf(x)`   | Elementwise log of standard normal CDF  |
| `cummax_expr(x)`   | Cumulative maximum along an axis        |
| `dotsort(X, W)`    | Weighted sorted dot product             |

#### Math function dispatch

Standard R math functions work directly on CVXR expressions:

`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``3``)`` `[`abs`](https://rdrr.io/r/base/MathFun.html)`(``x``)`` ``# elementwise absolute value`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``x``)`` ``# elementwise square root`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``x``)`` ``# sum of entries`` `[`max`](https://rdrr.io/r/base/Extremes.html)`(``x``)`` ``# maximum entry`` `[`norm`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md)`(``x``, ``"2"``)`` ``# Euclidean norm`

#### Boolean logic atoms

For mixed-integer programming:
[`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
[`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
[`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
[`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
[`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
[`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md).

#### Other new atoms

- `perspective(f, s)` for perspective functions
- `FiniteSet(expr, values)` constraint for discrete optimization
- [`ceil_expr()`](https://www.cvxgrp.org/CVXR/reference/ceil_expr.md),
  [`floor_expr()`](https://www.cvxgrp.org/CVXR/reference/floor_expr.md)
  for DQCP problems
- [`condition_number()`](https://www.cvxgrp.org/CVXR/reference/condition_number.md),
  [`gen_lambda_max()`](https://www.cvxgrp.org/CVXR/reference/gen_lambda_max.md),
  [`dist_ratio()`](https://www.cvxgrp.org/CVXR/reference/dist_ratio.md)
  for DQCP

### Backward-Compatibility Aliases

- [`tv()`](https://www.cvxgrp.org/CVXR/reference/tv.md) is
  **deprecated**; use
  [`total_variation()`](https://www.cvxgrp.org/CVXR/reference/total_variation.md)
  (still works but warns once)
- `norm2(x)` is **deprecated**; use `p_norm(x, 2)` (still works but
  warns once)
- `multiply(x, y)` is **deprecated**; use `x * y` for elementwise
  multiplication
- Old [`solve()`](https://rdrr.io/r/base/solve.html) still works and
  returns a compatibility list
- Old function names (`problem_status`, `getValue`, etc.) still work but
  emit once-per-session deprecation warnings

### Migration Guide

To migrate code from CVXR 1.x to 1.8.x:

1.  Replace `result <- solve(problem)` with `opt_val <- psolve(problem)`

2.  Replace `result$getValue(x)` with `value(x)`

3.  Replace `result$value` with the return value from
    [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md)

4.  Replace `result$status` with `status(problem)`

5.  Replace `result$getDualValue(con)` with `dual_value(con)`

6.  Update solver names: `"ECOS"` → `"CLARABEL"`, `"GLPK"` → `"HIGHS"`

7.  Update `axis` arguments: `axis = NA` → `axis = NULL` (row/column
    axis values 1 and 2 are unchanged)

8.  Replace `A %>>% B` with `PSD(A - B)` if desired

9.  Wrap Matrix package objects with
    [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md)
    before using them in CVXR expressions (e.g., `as_cvxr_expr(A) %*% x`
    instead of `A %*% x` when `A` is a `dgCMatrix` or other Matrix
    class). This preserves sparsity. Base R matrices need no wrapping.

10. **Dimension-preserving operations.** CVXR 1.8 preserves 2D shapes
    throughout, matching CVXPY. In particular, axis reductions like
    `sum_entries(X, axis = 2)` now return a proper row vector of shape
    `(1, n)` rather than collapsing to a 1D vector. When comparing such
    a result with an R numeric vector (which CVXR treats as a column),
    you may need to use [`t()`](https://rdrr.io/r/base/t.html) or
    `matrix(..., nrow = 1)` to match shapes:

    `## Old (worked in CVXR 1.x because axis reductions were 1D):`` `[`sum_entries`](https://www.cvxgrp.org/CVXR/reference/sum_entries.md)`(``X``, axis ``=`` ``2``)`` ``==`` ``target_vec`` ``## New (wrap target as row vector to match the (1, n) shape):`` `[`sum_entries`](https://www.cvxgrp.org/CVXR/reference/sum_entries.md)`(``X``, axis ``=`` ``2``)`` ``==`` `[`t`](https://rdrr.io/r/base/t.html)`(``target_vec``)`

    Similarly, if you extract a scalar from a CVXR result and need a
    plain numeric value, use
    [`as.numeric()`](https://rdrr.io/r/base/numeric.html) to drop the
    matrix dimensions.

### CRAN Submission Tip

If you encounter issues involving the `Rmosek` package while submitting
your package to CRAN, include the following code in `<your_pkg>/R/zzz.R`
to resolve the issue.

`## Content of <your_pkg>/R/zzz.R`` `` ``.onLoad`` ``<-`` ``function``(``libname``, ``pkgname``)`` ``{`` `` ``CVXR``::`[`exclude_solvers`](https://www.cvxgrp.org/CVXR/reference/available_solvers.md)`(``"MOSEK"``)`` ``}`` ``.onUnload`` ``<-`` ``function``(``libname``, ``pkgname``)`` ``{`` `` ``CVXR``::`[`include_solvers`](https://www.cvxgrp.org/CVXR/reference/available_solvers.md)`(``"MOSEK"``)`` ``}`

### Further Reading

- [CVXR website](https://cvxr.rbind.io) — worked examples
- [Package reference](https://www.cvxgrp.org/CVXR/) — full API
  documentation
- [CVXPY documentation](https://www.cvxpy.org/) — mathematical framework
- Fu, Narasimhan, and Boyd (2020). “CVXR: An R Package for Disciplined
  Convex Optimization.” *Journal of Statistical Software*, 94(14).
  <doi:10.18637/jss.v094.i14>
