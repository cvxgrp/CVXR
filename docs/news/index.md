# Changelog

## CVXR 1.9.1

This is the first CRAN release since 1.8.2 and is a large one: it tracks
CVXPY 1.9.1 and folds in the changes from the internal 1.8.2-1 and 1.9.0
development cycles. The headline additions are a derivative API for
differentiable convex programs, disciplined nonlinear programming
(DNLP), and interval-bounds propagation with native solver-bound
support.

### Disciplined nonlinear programming (DNLP)

- Solve smooth nonlinear programs with `psolve(prob, nlp = TRUE)`. The
  problem must satisfy the new
  [`is_dnlp()`](https://www.cvxgrp.org/CVXR/reference/is_dnlp.md)
  grammar (disciplined nonlinear program); DCP problems are a subset, so
  any DCP problem is also a valid DNLP.
- New smooth atoms usable in DNLPs:
  [`sin()`](https://rdrr.io/r/base/Trig.html),
  [`cos()`](https://rdrr.io/r/base/Trig.html),
  [`tan()`](https://rdrr.io/r/base/Trig.html),
  [`sinh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`asinh()`](https://rdrr.io/r/base/Hyperbolic.html),
  [`atanh()`](https://rdrr.io/r/base/Hyperbolic.html), and
  [`normcdf()`](https://www.cvxgrp.org/CVXR/reference/normcdf.md).
- [`prod()`](https://rdrr.io/r/base/prod.html) is now recognized as a
  smooth atom, so problems involving products of entries (including
  `prod_entries(X, axis = ...)`) are DNLP and solvable through the NLP
  path.
- The classically-differentiable atoms now report
  [`is_atom_smooth()`](https://www.cvxgrp.org/CVXR/reference/is_atom_smooth.md)
  matching CVXPY 1.9: affine atoms, `exp`, `log`, `entr`, `logistic`,
  `kl_div`, `rel_entr`, `xexp`, `power` (with a constant exponent),
  `geo_mean`, `log_sum_exp`, `quad_form`, `quad_over_lin`, and `prod`.
  As a result smooth convex/concave functions such as
  [`sum_squares()`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)
  are linearizable in both directions
  ([`is_smooth()`](https://www.cvxgrp.org/CVXR/reference/is_smooth.md)
  is `TRUE`), while piecewise atoms (`abs`, `min`, `max`) remain
  convex/concave-only. New predicate exports:
  [`is_smooth()`](https://www.cvxgrp.org/CVXR/reference/is_smooth.md),
  [`is_atom_smooth()`](https://www.cvxgrp.org/CVXR/reference/is_atom_smooth.md),
  [`is_linearizable_convex()`](https://www.cvxgrp.org/CVXR/reference/is_linearizable_convex.md),
  [`is_linearizable_concave()`](https://www.cvxgrp.org/CVXR/reference/is_linearizable_concave.md).
- Piecewise-linear atoms used in a DNLP (`abs`, `max`, `max_elemwise`,
  `norm_inf`, `sum_largest`, and `sum_smallest`) now initialize their
  epigraph variables during canonicalization, so the NLP backend has a
  complete starting point. `sum_largest` initializes its threshold to
  the (k+1)-th largest entry, matching CVXPY’s warm-start fix.
- New [`convolve()`](https://www.cvxgrp.org/CVXR/reference/convolve.md)
  atom, CVXPY 1.9’s preferred (numpy-style) name for the 1-D
  discrete-convolution
  [`conv()`](https://www.cvxgrp.org/CVXR/reference/conv.md) atom. For
  numeric input it falls through to
  [`stats::convolve()`](https://rdrr.io/r/stats/convolve.html).
- The DNLP path is powered by derivatives from the optional `sparsediff`
  package and can solve through the optional `ipopt` and `Uno` R
  packages. NLP solver names: `"IPOPT"`, `"UNO"` (and the variants
  `"uno_ipm"` / `"uno_sqp"`), `"KNITRO"`, and `"COPT"`; `"KNITRO"` and
  `"COPT"` are registered but require solver bindings not yet available
  in R. When installed, `IPOPT` is the first automatic NLP solver,
  matching CVXPY’s preference order.
- [`quad_over_lin()`](https://www.cvxgrp.org/CVXR/reference/quad_over_lin.md)
  and
  [`sum_squares()`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)
  now canonicalize axis-aware reductions (`axis = 1` and `axis = 2`)
  into batched second-order cone constraints, closing the CVXPY parity
  gap for `sum_squares(..., axis = ...)` in objectives and constraints.
- **The `UNO` interface defaults to the interior-point `ipopt` preset
  (plus MUMPS), which differs from CVXPY’s `filtersqp` default.**
  CVXPY’s `filtersqp` uses the BQPD active-set solver for its quadratic
  subproblems, which handles indefinite Hessians. The bundled Uno build
  is HiGHS-only (BQPD is not bundled and is not CRAN-license
  compatible), and HiGHS solves only *convex* quadratic subproblems, so
  `filtersqp` fails on nonconvex DNLPs. The interior-point `ipopt`
  preset factors the regularized KKT system with MUMPS and is robust on
  both convex and nonconvex problems. Force the SQP path with
  `solver = "uno_sqp"` (or `preset = "filtersqp"`); it is reliable only
  for problems whose subproblem Hessians stay positive semidefinite.
- Best-of-N random restarts for nonconvex DNLPs:
  `psolve(prob, nlp = TRUE, best_of = n)` solves from `n` random initial
  points and keeps the best result. Random initialization draws from a
  variable’s `sample_bounds(var) <- c(low, high)` (when set) or its
  finite variable bounds. The per-run objectives are available via
  `solver_stats(prob)@extra_stats$all_objs_from_best_of`.
- After an NLP solve,
  [`dual_value()`](https://www.cvxgrp.org/CVXR/reference/dual_value.md)
  returns the constraint duals recovered from the solver (a CVXR
  addition; the duals are not exposed by CVXPY’s NLP interface).

### Derivative API for differentiable convex programs

- New derivative API backed by the `diffcp` R package:
  `psolve(prob, requires_grad = TRUE)` followed by `backward(prob)` or
  `derivative(prob)`; `gradient(x)<-` and `delta(x)<-` set tangent /
  cotangent values on leaves. The chain rule is wired through `Dgp2Dcp`
  (log/exp) and `Complex2Real` (real/imag split).

### Bounds propagation

- [`get_bounds()`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)
  now works on any expression, not just variables: it propagates
  interval bounds through affine, elementwise, and piecewise-linear
  atoms (e.g. `get_bounds(A %*% x + b)`, `get_bounds(abs(x))`,
  `get_bounds(sum(x))`). Bounds are returned shaped to the expression’s
  dimensions.
- Variable bounds may now be sparse `Matrix` objects or symbolic bounds
  involving `Parameter`s.
  [`get_bounds()`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)
  preserves sparse numeric bounds and skips symbolic bounds, while
  solve-time attribute lowering enforces the symbolic constraints.
  `Parameter` values containing matching `Inf` entries now validate
  correctly.
- Solvers that expose native variable bounds now consume dense numeric
  bounds directly instead of adding bound constraints: HiGHS (LP/MILP
  and QP paths), Gurobi (QP/MIQP and conic/SOCP), CPLEX (QP/MIQP),
  XPRESS (QP/MIQP and conic/SOCP), PIQP (QP), and SCIP (conic). Sparse
  bounds continue through constraint lowering.
- Parametric variable bounds use bound tensors (HiGHS), so changing
  `Parameter` bounds updates native solver bounds on DPP re-solves.
- Bounds inferred for DCP auxiliary variables are now preserved for
  solvers that consume native variable bounds.
- Attribute reduction now compacts 2D symmetric, PSD, NSD, and diagonal
  `Parameter`s before rebuilding their full expressions.
- Variables now report `Parameter`s embedded in expression bounds and
  include those bounds in DPP/DGP compliance checks.

### Geometric and parameterized programming

- Positive (DGP) variables now accept numeric *and* parametric bounds
  under `gp = TRUE` (e.g. `Variable(pos = TRUE, bounds = list(lb, ub))`
  with `lb`/`ub` `Parameter`s). The DGP reduction log-transforms the
  bounds into the log domain and lowers them to constraints; parametric
  bounds canonicalize through the DGP tree without eagerly evaluating
  `log(value(param))`, so DPP re-solves with changed bound parameters
  work.
- [`is_dpp()`](https://www.cvxgrp.org/CVXR/reference/is_dpp.md) gains a
  `context` argument (`"dcp"` or `"dgp"`), matching CVXPY’s
  `is_dpp(context=...)`. In the `"dgp"` context a variable’s symbolic
  bounds must be log-log-affine (e.g. a product of positive `Parameter`s
  is DGP-DPP though not DCP-DPP;
  [`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md)/[`sum()`](https://rdrr.io/r/base/sum.html)
  bounds are rejected).
- New
  [`partial_optimize()`](https://www.cvxgrp.org/CVXR/reference/partial_optimize.md)
  transform.
- CPLEX now solves LP/SOCP/MI-LP/MI-SOCP problems through the conic path
  by translating SOC blocks to `Rcplex` quadratic constraints.

### New atoms, transforms, and accessors

- New [`sign()`](https://rdrr.io/r/base/sign.html) atom (DQCP, scalar
  input).
- [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) gains a
  `solver_path` argument for solver fallback chains.
- New
  [`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md) /
  `label<-` /
  [`format_labeled()`](https://www.cvxgrp.org/CVXR/reference/format_labeled.md)
  for attaching user-supplied labels to expressions.
- [`power()`](https://www.cvxgrp.org/CVXR/reference/power.md),
  [`geo_mean()`](https://www.cvxgrp.org/CVXR/reference/geo_mean.md), and
  [`p_norm()`](https://www.cvxgrp.org/CVXR/reference/p_norm.md) with
  `approx = TRUE` now warn when the selected solver supports power cones
  (`approx = FALSE` uses the exact `PowCone3D` / `PowConeND` form).
- MOSEK now supports mixed-integer programs (requires `Rmosek` 11.1.1+).
- [`solver_stats()`](https://www.cvxgrp.org/CVXR/reference/solver_stats.md)
  for MOSEK now populates `solve_time`, `num_iters`, and `extra_stats`.
- New
  [`param_dict()`](https://www.cvxgrp.org/CVXR/reference/param_dict.md),
  [`var_dict()`](https://www.cvxgrp.org/CVXR/reference/var_dict.md), and
  [`size_metrics()`](https://www.cvxgrp.org/CVXR/reference/size_metrics.md)
  accessors on `Problem`.
- `Variable(value = ...)` accepts an initial value at construction.
- New `CallbackParam` class — a `Parameter` subclass whose value is
  computed on every read by a user-supplied callback.
- [`project()`](https://www.cvxgrp.org/CVXR/reference/project.md)
  accepts index-based `boolean = c(i, j, ...)` and
  `integer = c(i, j, ...)` attributes (1-based).
- [`validate_arguments()`](https://www.cvxgrp.org/CVXR/reference/validate_arguments.md)
  now called on `Cummax`, `Cumprod`, and `MinEntries` constructors.

### Performance

- Problem canonicalization and solving are faster than the 1.8.2 CRAN
  release on a wide range of problems. Two changes account for most of
  the gain: a new `.fast_new` helper for S7 expression construction
  (roughly 40% lower wall-clock on atom-dense and PSD-heavy problems),
  and routing S7 class-membership checks on the hot canonicalization
  path through a cached check on the object’s class vector instead of
  [`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
  (which rebuilt the class name on every call; the per-node
  cone-detection scan benefits most). Deterministic memory allocation is
  unchanged. These are internal changes with no user-visible API
  difference.

### Bug fixes

- Bug fix (also affected 1.8.x):
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
  and
  [`get_problem_data()`](https://www.cvxgrp.org/CVXR/reference/get_problem_data.md)
  now take an explicit `gp` argument. Previously `gp = TRUE` passed
  through `...` was silently ignored, compiling a geometric program as a
  DCP problem.
- Bug fix (also affected 1.8.x):
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) now
  takes explicit `enforce_dpp` and `ignore_dpp` arguments, matching
  CVXPY’s [`solve()`](https://rdrr.io/r/base/solve.html). Previously
  they were silently swallowed by `...`; `enforce_dpp = TRUE` now raises
  on a non-DPP parametrized problem instead of falling back to a non-DPP
  compile.
- Fixed QP-path canonicalization of `quad_over_lin(expr, constant)` so
  the decision uses the canonicalized denominator, matching CVXPY 1.9
  behavior for cases such as `quad_over_lin(exp(x), 1)` that should keep
  `ExpCone` without adding an SOC block.
- Fix `Variable(c(n, n), diag = TRUE)` handling through
  `CvxAttr2Constr`.
- Fix
  [`perspective()`](https://www.cvxgrp.org/CVXR/reference/perspective.md)
  canonicalizer on PSD, NSD, and diag matrix variables.
- Fix [`project()`](https://www.cvxgrp.org/CVXR/reference/project.md) on
  sparse `Matrix`-package objects.

### Notes

- Complex-parameter DPP fast-path parity remains a deliberate feature
  gap: complex parameters are evaluated before `Complex2Real` because
  CVXR’s R sparse-matrix backend cannot carry complex sparse parameter
  tensors. Ordinary complex expression solving is still supported
  through `Complex2Real`.
- CVXR’s 2-D positive-dimensional model does not support zero-sized
  expressions, constraints, or solves. Generated or data-dependent
  models should skip vacuous constraints and use explicit zero
  contributions when a mask or filter selects no entries.

## CVXR 1.8.2

CRAN release: 2026-04-04

### New solvers: SCIP and XPRESS (15 total)

- Added SCIP solver (14th solver) via the R `scip` package. Supports LP,
  SOCP, MI-LP, and MI-SOCP problems. SCIP is registered in the conic
  solver path with `SUPPORTED_CONSTRAINTS = list(Zero, NonNeg, SOC)` and
  `MIP_CAPABLE = TRUE`.
- SCIP solver parameters can be passed via `scip_params` sub-list
  (matching CVXPY convention) for path-style parameter names like
  `"limits/time"`, `"limits/gap"`.
- Duals are not currently extracted (R `scip` package lacks dual API).
- 29 tests mirroring CVXPY’s `TestSCIP` class.
- Added XPRESS solver (15th solver) via the R `xpress` package with
  dual-interface architecture: QP path (`XPRESS_QP_Solver`) for LP/QP
  and conic path (`XPRESS_Conic_Solver`) for SOCP/MI-LP/MI-SOCP. Both
  paths support Zero and NonNeg constraints; the conic path also
  supports SOC. MIP warm-start is supported.

### `diag()` and `norm()` dispatch fixes

- [`diag()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) now
  works on CVXR expressions: `diag(vector_expr)` creates a diagonal
  matrix (`DiagVec`), `diag(square_matrix_expr)` extracts the diagonal
  (`DiagMat`). Falls through to `Matrix::diag()` for non-Expression
  inputs, correctly handling both Matrix S4 objects and base R matrices.
- [`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md)
  fallthrough now delegates to
  [`Matrix::norm()`](https://rdrr.io/pkg/Matrix/man/norm-methods.html)
  instead of [`base::norm()`](https://rdrr.io/r/base/norm.html), fixing
  `norm(sparse_matrix)` when CVXR is loaded.
- [`t()`](https://rdrr.io/r/base/t.html) and
  [`mean()`](https://rdrr.io/r/base/mean.html) switched from S7
  `method()` to
  [`registerS3method()`](https://rdrr.io/r/base/ns-internal.html) to
  avoid demoting Matrix’s S4 generics.

### CVXPY 1.8.2 parity

All applicable bug fixes from CVXPY v1.8.2 have been ported, each
annotated with `## CVXPY v1.8.2 fix:` in the source.

### `solver_opts()` and `use_quad_obj`

- New
  [`solver_opts()`](https://www.cvxgrp.org/CVXR/reference/solver_opts.md)
  constructor for unified solver options. Standard tolerance parameters
  (`feastol`, `reltol`, `abstol`, `num_iter`) and solver-specific
  parameters now flow through `psolve(...)` via the
  [`solver_opts()`](https://www.cvxgrp.org/CVXR/reference/solver_opts.md)
  constructor.
- New `use_quad_obj` option (default `TRUE`): when `FALSE`, forces conic
  decomposition path instead of QP, enabling `quad_form_canon` to detect
  indefinite P matrices.
- New
  [`supports_quad_obj()`](https://www.cvxgrp.org/CVXR/reference/supports_quad_obj.md)
  generic on conic solvers (Clarabel and SCS return `TRUE`, others
  `FALSE`).

### Element-wise matrix indexing (`SpecialIndex`)

- New `SpecialIndex` atom (mirroring CVXPY’s `special_index`) enables
  R-idiomatic element-wise selection on matrix expressions:

  - 2-column matrix indexing: `x[cbind(rows, cols)]`
  - Logical matrix indexing: `x[mask]`
  - Linear integer indexing: `x[c(1, 5, 9)]`
  - Logical vector indexing: `x[c(TRUE, FALSE, ...)]`

- This makes it natural to constrain specific entries of a matrix
  variable, e.g., fixing observed entries of a partially-known matrix:

  ``` r

  ind <- which(!is.na(Rmiss), arr.ind = TRUE)
  prob <- Problem(Minimize(obj), list(X[ind] == Rmiss[ind]))
  ```

- Previously, `x[i]` on a matrix expression errored with “Single-index
  selection on matrices not supported.” The workaround required an
  element-wise multiply with a binary mask, which was unintuitive.

- Internally uses sparse selection matrix multiplication (reshape →
  sparse matmul), requiring no C++ changes.

### Bug Fixes

- Fixed O(n²) memory usage in
  [`sum_squares()`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md),
  `power(x, 2)`,
  [`quad_over_lin()`](https://www.cvxgrp.org/CVXR/reference/quad_over_lin.md),
  and [`huber()`](https://www.cvxgrp.org/CVXR/reference/huber.md) when
  the first argument has many elements (e.g., regression residuals with
  thousands of observations). The quadratic canonicalizer was converting
  a sparse identity matrix to dense form. Memory now scales linearly
  with problem size.

## CVXR 1.8.1

CRAN release: 2026-03-06

### Complete rewrite using S7 object system

This is a ground-up rewrite of CVXR using R’s S7 object system, designed
to be isomorphic with CVXPY 1.8.1 for long-term maintainability. ~4-5x
faster than CVXR 1.0-15 on typical problems.

#### New features

- S7 class system replaces S4 for all expression, constraint, and
  problem classes. Significantly faster construction and method
  dispatch.
- 13 solvers at initial release: Clarabel (default), SCS, OSQP, HiGHS,
  MOSEK, Gurobi, GLPK, GLPK_MI, ECOS, ECOS_BB, CPLEX, CVXOPT, PIQP.
- Mixed-integer programming via GLPK_MI, ECOS_BB, Gurobi, or HiGHS
  (`boolean = TRUE` or `integer = TRUE` in
  [`Variable()`](https://www.cvxgrp.org/CVXR/reference/Variable.md)).
- Parameter support via
  [`Parameter()`](https://www.cvxgrp.org/CVXR/reference/Parameter.md)
  class with DPP (Disciplined Parameterized Programming) for efficient
  re-solves when only parameter values change.
- 50+ atom classes covering LP, QP, SOCP, SDP, exponential cone, and
  power cone problems.
- [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) as the
  primary solve interface, returning the optimal value directly.
- [`solve()`](https://rdrr.io/r/base/solve.html) backward-compatibility
  wrapper returning a `cvxr_result` list with `$value`, `$status`,
  `$solver`, `$getValue()`, and `$getDualValue()`.
- `verbose = TRUE` option in
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) for
  structured solve output with timing information.
- Standard solver parameters (`feastol`, `reltol`, `abstol`, `num_iter`)
  in [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) with
  automatic translation to solver-native names via
  [`solver_default_param()`](https://www.cvxgrp.org/CVXR/reference/solver_default_param.md).
  Solver-native parameters in `...` take priority.
- Automatic solver selection based on problem type.
- Warm-start support for 7 solvers (OSQP, SCS, Clarabel, Gurobi, MOSEK,
  PIQP; HiGHS blocked by R package limitation).
- Decomposed solve API:
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md),
  [`solve_via_data()`](https://www.cvxgrp.org/CVXR/reference/solve_via_data.md),
  [`problem_unpack_results()`](https://www.cvxgrp.org/CVXR/reference/problem_unpack_results.md)
  for compile-once/solve-many workflows.
- [`visualize()`](https://www.cvxgrp.org/CVXR/reference/visualize.md)
  for problem introspection: text display with Smith form annotations,
  interactive HTML output (D3 tree + KaTeX), on-demand DCP violation
  analysis for non-compliant problems, curvature coloring on constraint
  nodes, and matrix stuffing visualization (Stages 4-5: Standard Form +
  Solver Data) via the new `solver` parameter.
- Matrix package interoperability via
  [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md).
  Sparse Matrix objects (`dgCMatrix`, `dgeMatrix`, etc.) use S4 dispatch
  which preempts S7/S3, so direct arithmetic like `sparseA %*% x` fails.
  Wrapping with `as_cvxr_expr(sparseA) %*% x` converts to a CVXR
  Constant while preserving sparsity. Base R `matrix` and `numeric` work
  natively without wrapping.
- Full 2D broadcasting for `+` and `-` operators, matching CVXPY’s
  `Expression.broadcast()`. Expressions with compatible but different
  shapes (e.g., `(1, n) + (m, 1)` → `(m, n)`) now work correctly in
  constraints and objectives. Previously only scalar promotion was
  supported.
- Correct solver-cone routing matching CVXPY’s atom classifications.
  Atoms that canonicalize to PSD cones (`MatrixFrac`, `TrInv`,
  `SigmaMax`, `NormNuc`, `LambdaSumLargest`) are now correctly
  classified, ensuring solvers like ECOS that don’t support SDP are
  rejected with a clear error message instead of a cryptic dimension
  mismatch. Approximate variants (`PnormApprox`, `PowerApprox`,
  `GeoMeanApprox`) correctly route to SOC; exact variants route to their
  native cone types (`PowCone3D`, `PowConeND`).

#### Breaking changes from CVXR 1.x

- S7 classes replace S4. Use `Variable(n)` instead of
  `new("Variable", n)`.
- Requires R \>= 4.3.0 (for
  [`chooseOpsMethod()`](https://rdrr.io/r/base/chooseOpsMethod.html) and
  S7 `@` support).
- [`solve()`](https://rdrr.io/r/base/solve.html) now returns a
  `cvxr_result` S3 list, not an S4 object. Use
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) for a
  direct numeric return value.
- `getValue(x)` and `getDualValue(con)` are deprecated. Use `value(x)`
  and `dual_value(con)` after solving instead.
- [`problem_status()`](https://www.cvxgrp.org/CVXR/reference/problem_status.md),
  [`problem_solution()`](https://www.cvxgrp.org/CVXR/reference/problem_solution.md),
  and
  [`get_problem_data()`](https://www.cvxgrp.org/CVXR/reference/get_problem_data.md)
  are deprecated. Use
  [`status()`](https://www.cvxgrp.org/CVXR/reference/status.md),
  [`solution()`](https://www.cvxgrp.org/CVXR/reference/solution.md), and
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
  instead.
- The `axis` parameter uses 1-based indexing matching R’s
  [`apply()`](https://rdrr.io/r/base/apply.html) MARGIN convention:
  `axis = 1` for row-wise (reduce columns), `axis = 2` for column-wise
  (reduce rows), `axis = NULL` for all entries.
- Internal class structure completely changed (S7 instead of S4). Code
  that accessed slots directly will need updating.
- S4 class checks like `is(x, "Variable")` no longer work. Use
  `S7_inherits(x, Variable)`.

#### Convenience atoms and functions

- `ptp(x, axis, keepdims)` – peak-to-peak (range): `max(x) - min(x)`.
- `cvxr_mean(x, axis, keepdims)` – arithmetic mean along an axis.
- `cvxr_std(x, axis, keepdims, ddof)` – standard deviation.
- `cvxr_var(x, axis, keepdims, ddof)` – variance.
- `vdot(x, y)` – vector dot product (inner product).
- `scalar_product(x, y)` – alias for `vdot`.
- `cvxr_outer(x, y)` – outer product of two vectors.
- `inv_prod(x)` – reciprocal of product of entries.
- `loggamma(x)` – elementwise log of gamma function (piecewise linear
  approximation).
- `log_normcdf(x)` – elementwise log of standard normal CDF (quadratic
  approximation).
- `cummax_expr(x, axis)` – cumulative maximum along an axis.
- `dotsort(X, W)` – weighted sorted dot product (generalization of
  `sum_largest`/`sum_smallest`).
- [`diff_pos()`](https://www.cvxgrp.org/CVXR/reference/diff_pos.md),
  [`resolvent()`](https://www.cvxgrp.org/CVXR/reference/resolvent.md) –
  DGP convenience functions.

#### Clean API names

- `status(prob)` – returns the solution status (replaces
  [`problem_status()`](https://www.cvxgrp.org/CVXR/reference/problem_status.md)).
- `solution(prob)` – returns the raw Solution object (replaces
  [`problem_solution()`](https://www.cvxgrp.org/CVXR/reference/problem_solution.md)).
- `problem_data(prob, solver)` – returns solver-ready data (replaces
  [`get_problem_data()`](https://www.cvxgrp.org/CVXR/reference/get_problem_data.md)).
- Old names still work but emit once-per-session deprecation warnings.

#### Backward-compatibility aliases

- [`tv()`](https://www.cvxgrp.org/CVXR/reference/tv.md) is a deprecated
  alias for
  [`total_variation()`](https://www.cvxgrp.org/CVXR/reference/total_variation.md).
- `norm2(x)` is a deprecated alias for `p_norm(x, 2)`.
- `multiply(x, y)` is a deprecated alias for `x * y`.
- [`installed_solvers()`](https://www.cvxgrp.org/CVXR/reference/installed_solvers.md)
  lists available solver packages.

#### Advanced features

- DPP (Disciplined Parameterized Programming) for efficient parameter
  re-solve with compilation caching.
- DGP (Disciplined Geometric Programming) via `psolve(prob, gp = TRUE)`.
  New atoms:
  [`prod_entries()`](https://www.cvxgrp.org/CVXR/reference/prod_entries.md),
  [`cumprod()`](https://rdrr.io/r/base/cumsum.html),
  [`one_minus_pos()`](https://www.cvxgrp.org/CVXR/reference/one_minus_pos.md),
  [`eye_minus_inv()`](https://www.cvxgrp.org/CVXR/reference/eye_minus_inv.md),
  [`pf_eigenvalue()`](https://www.cvxgrp.org/CVXR/reference/pf_eigenvalue.md),
  [`gmatmul()`](https://www.cvxgrp.org/CVXR/reference/gmatmul.md).
- DQCP (Disciplined Quasiconvex Programming) via
  `psolve(prob, qcp = TRUE)`. New atoms:
  [`ceil_expr()`](https://www.cvxgrp.org/CVXR/reference/ceil_expr.md),
  [`floor_expr()`](https://www.cvxgrp.org/CVXR/reference/floor_expr.md),
  [`condition_number()`](https://www.cvxgrp.org/CVXR/reference/condition_number.md),
  [`gen_lambda_max()`](https://www.cvxgrp.org/CVXR/reference/gen_lambda_max.md),
  [`dist_ratio()`](https://www.cvxgrp.org/CVXR/reference/dist_ratio.md).
- Complex variable support via `Variable(n, complex = TRUE)` and
  `Complex2Real` reduction. Atoms: `Conj_()`, `Real_()`, `Imag_()`,
  [`expr_H()`](https://www.cvxgrp.org/CVXR/reference/expr_H.md) for
  conjugate transpose.
- `perspective(f, s)` atom for perspective functions.
- `FiniteSet(expr, values)` constraint for discrete optimization.
- Boolean logic atoms:
  [`Not()`](https://www.cvxgrp.org/CVXR/reference/Not.md),
  [`And()`](https://www.cvxgrp.org/CVXR/reference/And.md),
  [`Or()`](https://www.cvxgrp.org/CVXR/reference/Or.md),
  [`Xor()`](https://www.cvxgrp.org/CVXR/reference/Xor.md),
  [`implies()`](https://www.cvxgrp.org/CVXR/reference/implies.md),
  [`iff()`](https://www.cvxgrp.org/CVXR/reference/iff.md).
- `%>>%` and `%<<%` operators for PSD and NSD constraints.
- [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md)
  helper for wrapping R objects as CVXR constants. Required for Matrix
  package objects (`dgCMatrix`, `dgeMatrix`, etc.) which use S4 dispatch
  that preempts S7/S3. Preserves sparsity (unlike
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html)). Base R
  `matrix`/`numeric` work natively without wrapping.

#### Bug fixes

- Fixed [`kron()`](https://www.cvxgrp.org/CVXR/reference/kron.md) KRON_L
  coefficient matrix construction.
- Fixed `CvxAttr2Constr` index type and matrix variable handling.
- Fixed Clarabel PowConeND bridge for uniquely-named cones.
- Fixed C++ `diag` offset for super/sub-diagonals (`k != 0`).
- Fixed DQCP bisection: infinity check, stale bounds for Maximize,
  conditional high-endpoint seeding.
- Fixed MOSEK dual extraction (`suc - slc` reconstruction).
- Fixed OSQP/HiGHS inequality dual sign convention.

#### Known limitations

- Matrix package objects (`dgCMatrix`, `dgeMatrix`, `ddiMatrix`,
  `sparseVector`) cannot be used directly with CVXR operators due to S4
  dispatch preempting S7/S3. Wrap with
  [`as_cvxr_expr()`](https://www.cvxgrp.org/CVXR/reference/as_cvxr_expr.md)
  first. Base R `matrix`/`numeric` work natively. This is an R dispatch
  limitation requiring upstream changes in S7 and/or Matrix.
- Warm-start for HiGHS is implemented (works with `highs` 1.14) but the
  upstream R `highs` PR exposing `setSolution()` has not yet been
  actioned by the maintainers; lives on branch `highs-warm-start`.
- Complex SOC dual variables are not recovered (CVXPY itself does not
  implement this — `complex2real.py:56` declares
  `UNIMPLEMENTED_COMPLEX_DUALS = (SOC, OpRelEntrConeQuad)`). Complex SOC
  *primal* canonicalization works.
- Deferred to a future release: RelEntrConeQuad / OpRelEntrConeQuad,
  quantum atoms, CPLEX conic (SOC) path.

## CVXR 1.0-15

CRAN release: 2024-11-07

- Revert `clarabel` requirement and use enhance rather than import until
  really necessary ([Issue
  142](https://github.com/cvxgrp/CVXR/issues/142)).
- Update return codes for `user_limit` etc to be `infeasible_inaccurate`
  to match [`CVXPY` for
  gurobi](https://github.com/cvxpy/cvxpy/pull/1270)
- Add an S3 print method for result from
  [`solve()`](https://rdrr.io/r/base/solve.html).
- Move `upper_tri_to_full` to C++
- Drop use of R6 classes
- Address bug in copying `Power` object ([Issue
  145](https://github.com/cvxgrp/CVXR/issues/142)).

## CVXR 1.0-14

CRAN release: 2024-06-27

- Address further inefficiency in use of diagonal matrices in
  `qp2quad_form.R`
- Add initial interface to clarabel solver

## CVXR 1.0-13

CRAN release: 2024-06-01

- Address inefficient processing of cones for MOSEK (Issue 137 reported
  by aszekMosek)
- Fix `extract_quadratic_coeffs` to use sparse matrix and sweep in place
  for better memory use (reported by Marissa Reitsma)

## CVXR 1.0-12

CRAN release: 2024-02-01

- `Rmosek` to be removed from CRAN, so moved to drat repo
- Cleaned up some problematic Rd files shown in CRAN results

## CVXR 1.0-11

CRAN release: 2022-10-30

- Being more careful about coercing to `dgCMatrix` via
  `(as(as(<matrix>, "CsparseMatrix"), "generalMatrix"))`
- Modify all class inheritance checks to use
  [`inherits()`](https://rdrr.io/r/base/class.html)

## CVXR 1.0-10

CRAN release: 2021-11-10

- Now requiring the updated scs 3.0 as an import

## CVXR 1.0-9

CRAN release: 2021-01-19

- Now importing ECOSolveR version 0.5.4 and higher
- Added fixes for Matrix 1.3 update
- Somewhat better typesetting for Power class documentation (Issue
  [\#86](https://github.com/cvxgrp/CVXR/issues/86))

## CVXR 1.0-8 to 1.0-3

CRAN release: 2020-09-13

- Conforming to CRAN suggestions on non-CRAN packages, but no actual
  changes.

## CVXR 1.0-2

- Added exponential cone support for `MOSEK` and uncommented associated
  test.
- Added JSS publication DOI.

## CVXR 1.0-1

CRAN release: 2020-04-02

- Many small fixes to solver interfaces
- Reference semantics for `Parameter`
- Warm start for `OSQP` and updates to solver cache
- Solver parameter defaults explicit now
- New tests added for solver combinations

## CVXR 1.0

CRAN release: 2020-02-02

- Major release implementing reductions, many new solvers.

## CVXR 0.99

CRAN release: 2018-05-26

- Bug fix: duplicated integer and boolean indices.
- Bug fix: correct typo in constraint specification to GLPK.
- Added tutorial articles based on v0.99 to [CVXR
  website](https://cvxr.rbind.io) on using other solvers, integer
  programming, MOSEK and GUROBI examples.

## CVXR 0.98-1

- Minor typographical fixes.

## CVXR 0.98

- Dropped `delay_load` parameter dropped in
  [`reticulate::import_from_path`](https://rstudio.github.io/reticulate/reference/import.html),
  per changes in `reticulate`.

- Cleaned up hooks into reticulate for commercial solvers.

## CVXR 0.97-1

- Minor typo and documentation fixes.

## CVXR 0.97

- Added `LPSOLVE` via
  [`lpSolveAPI`](https://cran.r-project.org/package=lpSolveAPI)
- Added `GLPK` via [`Rglpk`](https://cran.r-project.org/package=Rglpk)
- Added `MOSEK`
- Added `GUROBI`
- Bug fix: [issue](https://github.com/cvxgrp/CVXR/issues/25)
  [\#25](https://github.com/cvxgrp/CVXR/issues/25). All CVXR expressions
  retain dimensions. Culprit was `drop = FALSE` (in function
  `Index.get_special_slice`) as suspected.

## CVXR 0.96

- Added a note that CVXR can probably be compiled from source for
  earlier versions of R. This is
  [issue](https://github.com/cvxgrp/CVXR/issues/24)
  [\#24](https://github.com/cvxgrp/CVXR/issues/24)

- Using [pkgdown](https://pkgdown.r-lib.org). This also addresses
  [issue](https://github.com/cvxgrp/CVXR/issues/23)
  [\#23](https://github.com/cvxgrp/CVXR/issues/23)

- Bug fix: [issue](https://github.com/cvxgrp/CVXR/issues/28)
  [\#28](https://github.com/cvxgrp/CVXR/issues/28) Function `intf_sign`
  (`interface.R`) was unnecessarily using a tolerance parameter, now
  eliminated.

## CVXR 0.95

CRAN release: 2018-02-20

- Updated Solver.solve to adapt to new ECOSolveR. Require version 0.4 of
  ECOSolveR now.

- Updated `unpack_results` to behave exactly like in CVXPY. Added
  documentation and testthat tests. Documented in the
  [Speed](https://cvxr.rbind.io/examples/solvers/speed).

## CVXR 0.94-4

CRAN release: 2017-11-20

- First CRAN release 2017-11-20.
