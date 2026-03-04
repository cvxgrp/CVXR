# CVXR 1.8.1

## Complete rewrite using S7 object system

This is a ground-up rewrite of CVXR using R's S7 object system,
designed to be isomorphic with CVXPY 1.8.1 for long-term
maintainability. ~4-5x faster than CVXR 1.0-15 on typical problems.

### New features

* S7 class system replaces S4 for all expression, constraint, and
  problem classes. Significantly faster construction and method
  dispatch.
* 13 solvers: Clarabel (default), SCS, OSQP, HiGHS, MOSEK, Gurobi,
  GLPK, GLPK_MI, ECOS, ECOS_BB, CPLEX, CVXOPT, PIQP.
* Mixed-integer programming via GLPK_MI, ECOS_BB, Gurobi, or HiGHS
  (`boolean = TRUE` or `integer = TRUE` in `Variable()`).
* Parameter support via `Parameter()` class with DPP (Disciplined
  Parameterized Programming) for efficient re-solves when only
  parameter values change.
* 50+ atom classes covering LP, QP, SOCP, SDP, exponential cone, and
  power cone problems.
* `psolve()` as the primary solve interface, returning the optimal
  value directly.
* `solve()` backward-compatibility wrapper returning a `cvxr_result`
  list with `$value`, `$status`, `$solver`, `$getValue()`, and
  `$getDualValue()`.
* `verbose = TRUE` option in `psolve()` for structured solve output
  with timing information.
* Standard solver parameters (`feastol`, `reltol`, `abstol`, `num_iter`)
  in `psolve()` with automatic translation to solver-native names via
  `solver_default_param()`. Solver-native parameters in `...` take priority.
* Automatic solver selection based on problem type.
* Warm-start support for 7 solvers (OSQP, SCS, Clarabel, Gurobi,
  MOSEK, PIQP; HiGHS blocked by R package limitation).
* Decomposed solve API: `problem_data()`, `solve_via_data()`,
  `problem_unpack_results()` for compile-once/solve-many workflows.
* `visualize()` for problem introspection: text display with Smith
  form annotations, interactive HTML output (D3 tree + KaTeX),
  on-demand DCP violation analysis for non-compliant problems,
  curvature coloring on constraint nodes, and matrix stuffing
  visualization (Stages 4-5: Standard Form + Solver Data) via the
  new `solver` parameter.
* Matrix package interoperability via `as_cvxr_expr()`. Sparse Matrix
  objects (`dgCMatrix`, `dgeMatrix`, etc.) use S4 dispatch which
  preempts S7/S3, so direct arithmetic like `sparseA %*% x` fails.
  Wrapping with `as_cvxr_expr(sparseA) %*% x` converts to a CVXR
  Constant while preserving sparsity. Base R `matrix` and `numeric`
  work natively without wrapping.

### Breaking changes from CVXR 1.x

* S7 classes replace S4. Use `Variable(n)` instead of
  `new("Variable", n)`.
* Requires R >= 4.3.0 (for `chooseOpsMethod()` and S7 `@` support).
* `solve()` now returns a `cvxr_result` S3 list, not an S4 object.
  Use `psolve()` for a direct numeric return value.
* `getValue(x)` and `getDualValue(con)` are deprecated. Use `value(x)`
  and `dual_value(con)` after solving instead.
* `problem_status()`, `problem_solution()`, and `get_problem_data()` are
  deprecated. Use `status()`, `solution()`, and `problem_data()` instead.
* The `axis` parameter uses 1-based indexing matching R's `apply()`
  MARGIN convention: `axis = 1` for row-wise (reduce columns),
  `axis = 2` for column-wise (reduce rows), `axis = NULL` for all
  entries.
* Internal class structure completely changed (S7 instead of S4).
  Code that accessed slots directly will need updating.
* S4 class checks like `is(x, "Variable")` no longer work. Use
  `S7_inherits(x, Variable)`.

### Convenience atoms and functions

* `ptp(x, axis, keepdims)` -- peak-to-peak (range): `max(x) - min(x)`.
* `cvxr_mean(x, axis, keepdims)` -- arithmetic mean along an axis.
* `cvxr_std(x, axis, keepdims, ddof)` -- standard deviation.
* `cvxr_var(x, axis, keepdims, ddof)` -- variance.
* `vdot(x, y)` -- vector dot product (inner product).
* `scalar_product(x, y)` -- alias for `vdot`.
* `cvxr_outer(x, y)` -- outer product of two vectors.
* `inv_prod(x)` -- reciprocal of product of entries.
* `loggamma(x)` -- elementwise log of gamma function (piecewise linear
  approximation).
* `log_normcdf(x)` -- elementwise log of standard normal CDF (quadratic
  approximation).
* `cummax_expr(x, axis)` -- cumulative maximum along an axis.
* `dotsort(X, W)` -- weighted sorted dot product (generalization of
  `sum_largest`/`sum_smallest`).
* `diff_pos()`, `resolvent()` -- DGP convenience functions.

### Clean API names

* `status(prob)` -- returns the solution status (replaces `problem_status()`).
* `solution(prob)` -- returns the raw Solution object (replaces `problem_solution()`).
* `problem_data(prob, solver)` -- returns solver-ready data (replaces `get_problem_data()`).
* Old names still work but emit once-per-session deprecation warnings.

### Backward-compatibility aliases

* `tv()` is a deprecated alias for `total_variation()`.
* `norm2(x)` is a deprecated alias for `p_norm(x, 2)`.
* `multiply(x, y)` is a deprecated alias for `x * y`.
* `installed_solvers()` lists available solver packages.

### Advanced features

* DPP (Disciplined Parameterized Programming) for efficient parameter
  re-solve with compilation caching.
* DGP (Disciplined Geometric Programming) via `psolve(prob, gp = TRUE)`.
  New atoms: `prod_entries()`, `cumprod()`, `one_minus_pos()`,
  `eye_minus_inv()`, `pf_eigenvalue()`, `gmatmul()`.
* DQCP (Disciplined Quasiconvex Programming) via `psolve(prob, qcp = TRUE)`.
  New atoms: `ceil_expr()`, `floor_expr()`, `condition_number()`,
  `gen_lambda_max()`, `dist_ratio()`.
* Complex variable support via `Variable(n, complex = TRUE)` and
  `Complex2Real` reduction. Atoms: `Conj_()`, `Real_()`, `Imag_()`,
  `expr_H()` for conjugate transpose.
* `perspective(f, s)` atom for perspective functions.
* `FiniteSet(expr, values)` constraint for discrete optimization.
* Boolean logic atoms: `Not()`, `And()`, `Or()`, `Xor()`, `implies()`, `iff()`.
* `%>>%` and `%<<%` operators for PSD and NSD constraints.
* `as_cvxr_expr()` helper for wrapping R objects as CVXR constants.
  Required for Matrix package objects (`dgCMatrix`, `dgeMatrix`, etc.)
  which use S4 dispatch that preempts S7/S3. Preserves sparsity
  (unlike `as.matrix()`). Base R `matrix`/`numeric` work natively
  without wrapping.

### Bug fixes

* Fixed `kron()` KRON_L coefficient matrix construction.
* Fixed `CvxAttr2Constr` index type and matrix variable handling.
* Fixed Clarabel PowConeND bridge for uniquely-named cones.
* Fixed C++ `diag` offset for super/sub-diagonals (`k != 0`).
* Fixed DQCP bisection: infinity check, stale bounds for Maximize,
  conditional high-endpoint seeding.
* Fixed MOSEK dual extraction (`suc - slc` reconstruction).
* Fixed OSQP/HiGHS inequality dual sign convention.

### Known limitations

* Matrix package objects (`dgCMatrix`, `dgeMatrix`, `ddiMatrix`,
  `sparseVector`) cannot be used directly with CVXR operators due to
  S4 dispatch preempting S7/S3. Wrap with `as_cvxr_expr()` first.
  Base R `matrix`/`numeric` work natively. This is an R dispatch
  limitation requiring upstream changes in S7 and/or Matrix.
* Warm-start for HiGHS is blocked (R `highs` package lacks `setSolution()` API).
* Derivative/sensitivity API deferred (requires `diffcp`, no R equivalent).
* Deferred to a future release: RelEntrConeQuad, quantum atoms, MOSEK MIP.

# CVXR 1.0-15

* Revert `clarabel` requirement and use enhance rather than import
  until really necessary ([Issue
  142](https://github.com/cvxgrp/CVXR/issues/142)).
* Update return codes for `user_limit` etc to be
  `infeasible_inaccurate` to match [`CVXPY` for gurobi](https://github.com/cvxpy/cvxpy/pull/1270)
* Add an S3 print method for result from `solve()`.
* Move `upper_tri_to_full` to C++
* Drop use of R6 classes
* Address bug in copying `Power` object ([Issue 145](https://github.com/cvxgrp/CVXR/issues/142)).

# CVXR 1.0-14

* Address further inefficiency in use of diagonal matrices in
  `qp2quad_form.R`
* Add initial interface to clarabel solver

# CVXR 1.0-13

* Address inefficient processing of cones for MOSEK (Issue 137
  reported by aszekMosek)
* Fix `extract_quadratic_coeffs` to use sparse matrix and sweep in
  place for better memory use (reported by Marissa Reitsma)

# CVXR 1.0-12

* `Rmosek` to be removed from CRAN, so moved to drat repo
* Cleaned up some problematic Rd files shown in CRAN results

# CVXR 1.0-11

* Being more careful about coercing to `dgCMatrix`
  via `(as(as(<matrix>, "CsparseMatrix"), "generalMatrix"))`
* Modify all class inheritance checks to use `inherits()`

# CVXR 1.0-10

* Now requiring the updated scs 3.0 as an import

# CVXR 1.0-9

* Now importing ECOSolveR version 0.5.4 and higher
* Added fixes for Matrix 1.3 update
* Somewhat better typesetting for Power class documentation (Issue #86)

# CVXR 1.0-8 to 1.0-3

* Conforming to CRAN suggestions on non-CRAN packages, but no actual
  changes.

# CVXR 1.0-2

* Added exponential cone support for `MOSEK` and uncommented
  associated test.
* Added JSS publication DOI.

# CVXR 1.0-1

* Many small fixes to solver interfaces
* Reference semantics for `Parameter`
* Warm start for `OSQP` and updates to solver cache
* Solver parameter defaults explicit now
* New tests added for solver combinations

# CVXR 1.0

* Major release implementing reductions, many new solvers.

# CVXR 0.99

* Bug fix: duplicated integer and boolean indices.
* Bug fix: correct typo in constraint specification to GLPK.
* Added tutorial articles based on v0.99 to [CVXR
website](https://cvxr.rbind.io) on using other solvers, integer
programming, MOSEK and GUROBI examples.

# CVXR 0.98-1

* Minor typographical fixes.

# CVXR 0.98

* Dropped `delay_load` parameter dropped in
  `reticulate::import_from_path`, per changes in `reticulate`.

* Cleaned up hooks into reticulate for commercial solvers.

# CVXR 0.97-1

* Minor typo and documentation fixes.

# CVXR 0.97

* Added `LPSOLVE` via [`lpSolveAPI`](https://cran.r-project.org/package=lpSolveAPI)
* Added `GLPK` via [`Rglpk`](https://cran.r-project.org/package=Rglpk)
* Added `MOSEK`
* Added `GUROBI`
* Bug fix: [issue #25](https://github.com/cvxgrp/CVXR/issues/25).
  All CVXR expressions retain dimensions. Culprit was `drop =
  FALSE` (in function `Index.get_special_slice`) as suspected.

# CVXR 0.96

* Added a note that CVXR can probably be compiled from source for
  earlier versions of R. This is [issue
  #24](https://github.com/cvxgrp/CVXR/issues/24)

* Using [pkgdown](https://pkgdown.r-lib.org). This also addresses
  [issue #23](https://github.com/cvxgrp/CVXR/issues/23)

* Bug fix: [issue #28](https://github.com/cvxgrp/CVXR/issues/28)
  Function `intf_sign` (`interface.R`) was unnecessarily using a
  tolerance parameter, now eliminated.

# CVXR 0.95

* Updated Solver.solve to adapt to new ECOSolveR. Require version 0.4 of
  ECOSolveR now.

* Updated `unpack_results` to behave exactly like in CVXPY. Added
  documentation and testthat tests. Documented in the [Speed](https://cvxr.rbind.io/examples/solvers/speed).

# CVXR 0.94-4

* First CRAN release 2017-11-20.
