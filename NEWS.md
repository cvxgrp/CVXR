# CVXR 1.9.3

## Changed results

These fix wrong answers, so code that relied on the old results will see
different numbers.

* GUROBI solved quadratic objectives with cross terms (e.g.
  `quad_form(w, Sigma)`) with every cross term halved. Gurobi now receives
  the full matrix, as in CVXPY. Present in 1.9.2.
* CVXOPT returned wrong answers on problems with no inequality constraints
  (feasible equality-constrained problems were reported infeasible) and
  reported unbounded problems as optimal. These now solve correctly or fail
  with a solver error, as in CVXPY.
* GUROBI second-order-cone dual values had the wrong sign; they now agree
  with CVXPY and with every other solver.
* `cvxr_norm(X, p)` and `norm(X, type)` now agree with base R's `norm()` and
  CVXPY on matrices: with no `axis`, `p = 1`, `2` and `Inf` are the maximum
  absolute column sum, the largest singular value and the maximum absolute
  row sum. They were entrywise; use `cvxr_norm(X, "fro")` or `norm1(X)` for
  those.
  - A row vector is a 1 x n matrix, as in base R, so its 1- and Inf-norms
    are the other way round from before.
  - Any other `p` on a matrix or row vector is now an error, "Unsupported
    matrix norm.", as in CVXPY; use `p_norm(x, p)` for the entrywise norm.
  - `norm()` on an expression accepts every base R type (`"O"`, `"M"` and
    `"E"` failed, and `"I"` was the largest absolute entry).
* `gmatmul(A, X)` with a `Parameter` in both `A` and `X` gave a wrong answer
  under `gp = TRUE`, and parameter substitution rebuilt it with `A` and `X`
  swapped. Both fixed (CVXPY 1.9.3 still has the swap).
* The cached DGP, DQCP and problem-level curvature queries ignored whether a
  DPP check was running, so a DGP problem with parameters could be reported
  DPP and solved with wrong solver data. All caches are now keyed on the DPP
  scope, as CVXPY's are.
* `perspective()` of a function of a `symmetric`, `PSD` or `NSD` matrix
  variable of size 3 or more returned a wrong optimum.
* `sum_entries(x, axis = 2)` of a one-row expression, and any product of a
  1 x 1 constant with a multi-column expression, gave wrong results or
  crashed R. Present in 1.9.2.
* `sum()`, `max()`, `min()` and `prod()` of an expression now honor `axis`
  and `keepdims`, as in CVXPY. A named `axis` was taken as one more operand,
  so `sum(x, axis = 2)` returned `sum(x + 2)`; other named arguments are
  now an error.
* `matrix_trace(A * B)` returned the trace of the matrix product (CVXPY
  1.9.3, #3492).
* `cumsum(X)` of a matrix is now a column of `length(X)` entries in
  column-major order, as base R's `cumsum()` returns; `value()` and the
  solver used to disagree.
* `lambda_sum_largest(X, k)` accepts a fractional `k`, as CVXPY does; `k`
  was silently truncated.
* Bounds through `min_elemwise()` and `max_elemwise()` with operands of
  different shapes now broadcast as NumPy does; a row operand was permuted.
* `derivative()` and `backward()` through a Hermitian variable or parameter
  of size 4 or more had imaginary entries transposed.
* `backward()` through a `sparsity` variable returned garbage, and on the NLP
  path the bounds of a `sparsity` variable were dropped. A recovered
  variable value is now projected onto its attribute (a `nonneg` -1e-9 reads
  back as 0), as in CVXPY.
* `psolve(solver = "CLARABEL", warm_start = TRUE)` re-solves from scratch
  when the cones changed (e.g. a parametric power-cone exponent), instead of
  returning the previous optimum. CVXPY 1.9.3 still has this defect.
* A complex variable declared `PSD`/`NSD` lost its cone (#3342), and one
  declared `diag` or with a `sparsity` pattern solved as dense (CVXPY 1.9.3
  still ignores those). Both are now enforced.

## New features

* `geo_mean()` accepts an argument of any shape and gains `axis` and
  `keepdims` (CVXPY 1.9.3, #3493).
* `cvxr_var()` supports `axis` and `keepdims`, like `cvxr_mean()` and
  `cvxr_std()` (#3496).
* `cvxr_norm(X, "nuc")` returns the nuclear norm (it errored).
* Exact `geo_mean(x, approx = FALSE)`, and any `PowConeND` constraint, now
  solve on SCS and MOSEK, not just Clarabel. `PowConeND()` with a row-vector
  `z`, or a 1 x n `W` with `axis = 1`, now works and has correct duals.
* `perspective(f, s)` where `f` needs a positive semidefinite cone (for
  example `lambda_max()` or `sigma_max()`) now solves on Clarabel, SCS and
  MOSEK.
* `grad()` works for `partial_optimize()` expressions and
  `lambda_sum_largest()`, and is correct for complex arguments of
  `norm_nuc()`, `sigma_max()`, `quad_form()`, `quad_over_lin()` and
  `matrix_frac()`.
* `indicator()` and `partial_optimize()` expressions now have bounds, as in
  CVXPY: `(0, Inf)` for an indicator and unbounded for a partial problem.
* A quadratic objective over a `symmetric`, `PSD`, `NSD`, `diag` or
  `sparsity` variable now solves on QP-capable solvers (#3384), and
  `quad_over_lin()` / `sum_squares()` with an `axis` take the QP path, so
  OSQP accepts them.
* DNLP (`psolve(..., nlp = TRUE)`) now needs `sparsediff (>= 0.6.1)`
  and consequences thereof. It also accepts `NonNeg()`, `NonPos()` and
  `Zero()` constraints, which it used to abort on.
* New `curvatures()`, `log_log_curvature()` and `is_diag()`, as CVXPY's
  `Expression.curvatures`, `Expression.log_log_curvature` and
  `Leaf.is_diag()`.
* `curvature()` now reports log-log and quasi curvature (e.g.
  `"LOG-LOG CONVEX"`, `"QUASILINEAR"`) instead of `"UNKNOWN"`, as in CVXPY.

## DPP, DGP and DQCP

* `psolve()` now checks that every `Parameter` has a value before it starts
  any solve, as CVXPY does.
* `psolve(gp = TRUE)` now checks variable bounds when deciding whether a
  parameterized problem is DPP, as CVXPY does, and requires a conic solver.
* `is_dpp()` is now one DCP pass inside the DPP scope, as in CVXPY, instead of
  that pass plus a recursion over every argument.
* These are no longer DPP, as in CVXPY: `quad_over_lin(x, p)` with a
  `Parameter` denominator (#2433), a `PowCone3D()`/`PowConeND()` with a
  `Parameter` exponent, and `dotsort(x, W)` with a non-affine parametric
  `W`. `gmatmul(A, X)` requires a parametric `A` to be the `Parameter`
  itself. `cumsum()` of a `Parameter` is DPP (#3398).
* `psolve(requires_grad = TRUE)` checks match CVXPY: a non-DPP problem is a
  `DPPError`, `qcp = TRUE` is refused, and `solver = "SCS"` routes to
  DIFFCP.
* DQCP bisection now follows CVXPY: no extra solve at the upper end of an
  automatically found interval, sub-solve warnings reach the user, a failed
  query point warns and bisection continues, and its aborts carry class
  `SolverError`.

## Solvers

* Solver selection now follows CVXPY 1.9: the solver is chosen from the
  problem's structure first, and its interface (QP or conic) follows. The
  commercial-first default covers QPs and mixed-integer QPs, SCIP solves
  mixed-integer QPs, a named HiGHS, Gurobi, CPLEX or Xpress solves an LP
  through its conic interface, and falling back to a non-default solver
  warns. ECOS_BB is never chosen unless named.
* When no installed (or non-excluded) solver can be used, the error now reads
  "No installed solver could handle this problem.", as in CVXPY.
* `psolve(solver_path = ...)` forwards every `psolve()` argument to each
  attempt (`enforce_dpp`, `ignore_dpp` and `nlp` were dropped), and only a
  `SolverError` moves on to the next solver, as in CVXPY. The conditions
  CVXPY raises as `SolverError` now carry that class in R too.
* MOSEK: the status comes from the solution status, as in CVXPY, so a
  termination code such as an iteration or time limit is no longer an error
  by itself. The simplex optimizers work, and `bfs = TRUE` and
  `save_file = <path>` are accepted.
* HiGHS: options may be nested in a `highs_options` list, as in CVXPY, and
  an option HiGHS rejects is an error instead of being dropped. A stop on
  `mip_max_improving_sols` returns the incumbent with status `user_limit`.
* Gurobi: `warm_start = TRUE` seeds the start point from the variables'
  values, `save_file` writes the model, and `env` is rejected, as in CVXPY.
  `solver_stats(prob)@num_iters` is 0 instead of `NULL` when presolve
  finished the problem.
* CPLEX: `reoptimize = TRUE` re-solves with presolve off on "infeasible or
  unbounded", as in CVXPY; `solve_time` is set; unknown options abort. A
  continuous problem with a second-order cone now returns dual values (all
  were missing), as in CVXPY.
* SCS: an infeasible solve exposes the infeasibility certificate through
  `dual_value()` (#3197), and `extra_stats` is the full SCS result.
* `solver_stats(prob)@extra_stats` now carries the raw solver result for
  OSQP, ECOS, ECOS_BB, Gurobi, DIFFCP, HiGHS and PIQP, as CVXPY does, and
  `solution(prob)` holds only the problem's own variables.
* PIQP: `backend = "dense"` builds a dense model (it was silently sparse).

## Errors, warnings and other changes

* `power(x, p)` now takes its sign, curvature, domain, value and printed name
  from `p` as written, as CVXPY does, rather than from the rational
  approximation used to solve it.
* Calling `NonPos(expr)` explicitly now warns, once per session, that it is
  deprecated in favor of operator overloading (`expr <= 0`) or
  `NonNeg(-expr)`, in CVXPY's words and with class
  `CvxpyDeprecationWarning`.
* Now errors, as in CVXPY: a negative `axis`; an `axis` with
  `cvxr_norm(X, "fro")` or `"nuc"` (#3501); a complex variable declared
  `nonneg`, `nonpos`, `pos`, `neg` or with `bounds`, which were silently
  ignored (#3342); `partial_optimize()` of a complex inner problem.
* Whether a constant matrix is symmetric, Hermitian or skew-symmetric is now
  decided with NumPy's `allclose` tolerance (relative 1e-5, absolute 1e-8), as
  in CVXPY; `quad_form(x, P)` used to reject a `P` that CVXPY accepts.
* Smaller alignments with CVXPY 1.9.3 in `get_bounds()`, `grad()` and
  `domain()` of several atoms, `dual_cone()`, `is_symmetric()` of `Re(X)`,
  `constr_size()`, `is_dpp()` with an unsupported context, interval
  bounds of matrix products, and `is_quadratic()` of `x / y` with a
  non-constant `y` (now `FALSE`). `value(constraint, tolerance = )` exposes the
  violation tolerance.
* `visualize()` no longer errors on a problem containing `sqrt()`,
  `inv_pos()` or a fractional `power()` (broken in 1.9.2).
* Problems with no constraint rows no longer trigger the Matrix deprecation
  warning `as(<ngCMatrix>, "dgCMatrix") is deprecated`, and Gurobi QPs no
  longer trigger the `dtCMatrix` one.
* Internal: a conic solver's `reduction_accepts()` now refuses a
  mixed-integer problem when the solver is not MIP-capable, as CVXPY's
  `ConicSolver.accepts` does. `psolve()` already refused such a solve while
  choosing the solver, so no solve changes. `constr_expr()` returns
  `lhs - rhs` for `==`, `<=` and `>=`, as CVXPY's `.expr` does.

## Performance

* The first solve of a parametric (DPP) problem compiles about twice as fast
  and allocates a small fraction of the memory it did: at 500 variables with
  CLARABEL, `sum_squares(w - q)` goes from 21 ms / 24 MB to 11 ms / 0.9 MB, an
  SOC problem from 17.5 ms / 12.7 MB to 12 ms / 0.9 MB, an LP with a 100 x 500
  dense matrix from 55 ms / 23 MB to 50 ms / 8 MB.
* A DPP problem whose objective is quadratic in the variables and depends on
  parameters (`sum_squares(y - q)`, `quad_form(y, P) + t(q) %*% y`, ...) now
  keeps its parameter program across solves, as in CVXPY.
* The compiled-chain cache is keyed as CVXPY's is, on the solver, `gp`,
  `ignore_dpp` and `use_quad_obj`; `enforce_dpp` no longer takes part.
* DQCP bisections are 7 percent and medium LPs 4 percent faster than in
  CVXR 1.9.2.
* `psolve(nlp = TRUE)` hands a dense constant matrix with at least 5 percent
  nonzeros to the derivative engine as dense, as CVXPY does: 1.7 to 2.6 times
  faster for a dense 1000 x 1000 `A`. Pass a sparse `Matrix` to keep the
  sparse form for a structured matrix.
* `psolve(nlp = TRUE)` no longer makes dense copies of a sparse constant
  matrix while building the problem: a solve with a 3000 x 3000 one allocated
  208 MB and now allocates 2 MB.

# CVXR 1.9.2

This release tracks CVXPY 1.9.2. It is chiefly a correctness release: a number
of bugs that returned plausible wrong answers without any error are fixed, and
dual variables now come back with the shape of their constraint.

## Incorrect results fixed

* A constant weight on `sum(square(x))` was silently dropped, so penalized
  least-squares (ridge) problems were solved with weight 1 (since 1.8.2).
  `sum_squares()`, `quad_form()`, and unweighted terms were unaffected.
* A failed solve poisoned the compiled-chain cache: repeating the identical
  call could silently return an answer computed with the previous options, or
  die with an unrelated `NULL` error instead of the real message.
* The `sparsity` variable attribute was accepted and then ignored — the
  variable solved as fully dense. It is now enforced; off-pattern entries come
  back as structural zeros, matching CVXPY.
* `gradient(p)` from the derivative API was exactly 0 whenever
  `value(p) == 0`, and correct elsewhere. Fixed jointly with `diffcp`;
  requires `diffcp (>= 0.1.2)`.
* DQCP problems with nested inverses could return a confidently wrong
  "optimal" answer: the bisection parameter was canonicalized away, so every
  query point looked feasible.
* `Variable(pos = TRUE)` and `neg = TRUE` emitted no sign constraint, so
  bounded problems were reported unbounded on most solvers (since 1.8.2).
* Partial `boolean`/`integer` attributes were silently ignored — the problem
  solved as its continuous relaxation and reported `"optimal"`. They are now
  honored and require a mixed-integer solver.
* Numeric bounds of the wrong shape on a symmetric leaf were silently
  truncated (HiGHS returned a feasible but suboptimal answer). Bounds are now
  shape-checked, and conflicting dimension-reducing attributes
  (`symmetric` + `diag`, `PSD` + `NSD`, ...) are rejected, as in CVXPY.
* `project()` and value assignment disagreed with the solver path on partial
  `boolean`/`integer` indices: one index spelling rejected valid values,
  another silently accepted invalid ones.
* A complex `quad_form()` with a negative semidefinite matrix passed
  `is_dcp()` and then failed to compile; it now solves. An indefinite matrix
  is no longer mislabeled PSD.
* Quasiconvex (`qcp = TRUE`) problems treated a matrix product `A %*% x` as
  an elementwise one; it is now rejected up front with CVXPY's message.
  Elementwise and scalar products are unaffected.
* `vdot()` and `scalar_product()` omitted the conjugation for complex
  arguments. Matching upstream, `vdot()` of a non-affine real argument is now
  of unknown curvature.
* `conv(expr, kernel)` — the constant kernel second — read monotonicity off
  the wrong argument and accepted a concave composition as convex.
* `sqrt()` and `inv_pos()` now canonicalize to second-order cones as in
  CVXPY, not power cones: `ECOS` and mixed-integer SOCP solvers accept them,
  and an `SCS` quasiconvex failure mode (`NaN` objective) is gone. Use
  `power(x, 0.5, approx = FALSE)` for the power-cone form.
* Duals of a PSD constraint on a Hermitian matrix were half their correct
  value; PSD duals were transposed for sizes 3 and up on
  `CLARABEL`/`SCS`/`MOSEK`; `PowCone3D` duals were transposed across cones in
  multi-cone constraints. Objectives and primal solutions were never affected.
* A complex `Parameter` never triggered the complex-to-real reduction, so
  problems CVXPY solves errored with "Inequality constraints cannot be
  complex."
* Gurobi: `SUBOPTIMAL` mapped to `optimal_inaccurate`, claiming a near-optimal
  solution that need not exist, and the QP path promoted `SOLVER_ERROR` to
  `OPTIMAL_INACCURATE`. Both now match CVXPY.
* CPLEX: the status table's labels — and the mappings that followed them —
  were shifted by one across the abort-status block. Statuses are now derived
  from `(solstat, pfeas, dfeas)` as upstream does.
* Smaller fixes: `floor()` could be off by one within solver tolerance of an
  integer; `log1p()` reported the bounds of `log(x)`; `get_bounds()` errored
  on fractional powers; a `hermitian` leaf could not hold a complex value;
  assigning `NaN` to a leaf raised an R internal error instead of the real
  message.

## Dual variables have their constraint's shape

* `dual_value()` returns a matrix shaped like its constraint — n-by-n for a
  PSD constraint, n-by-1 for a vector constraint — instead of a flat vector.
  Scalar constraints, `SOC`, and `ExpCone` are exempt, matching CVXPY.
* Linear indexing (`d[3]`) still works; code calling `length()` or `dim()` on
  a dual may need adjusting.

## Mixed-integer variable attributes

* `boolean` and `integer` accept four spellings: whole-variable (`TRUE`),
  flat 1-based positions, `cbind(row, col)` pairs, and a logical mask of the
  variable's shape.
* Index lists are validated at construction — out-of-range, fractional,
  duplicated, or `NA` indices raise. A partial list combined with
  `symmetric`/`PSD`/`NSD`/`diag` is rejected, as in CVXPY; logic atoms
  (`And`, `Or`, ...) accept partially-boolean variables.

## Other improvements

* Common-subexpression elimination: a subtree appearing more than once is
  canonicalized once. Transparent — solutions, objective values, and duals
  are unchanged, and only structurally identical subtrees are merged.
* `Minimize()`/`Maximize()` support `+`, `-`, `*`, `/` as in CVXPY; negation
  or a negative scale flips the sense, and mixing senses in a sum raises.
  `is_dpp()` rejects an unrecognized context instead of silently answering
  the `"dcp"` question.
* DQCP bisection compiles its one-parameter subproblem family once when the
  reduction emits no lazy constraints, instead of recanonicalizing at every
  query — about 14x on `maximize sqrt(x)/exp(x)`. Answers are unchanged.
  An explicit bisection bound of 0 can now grow, and two valid quasiconvex
  forms previously rejected now solve.
* `norm(x, 1)` is recognized as log-log convex; `partial_trace()` and
  `partial_transpose()` propagate PSD/symmetric/hermitian structure;
  `atanh()` reports its domain instead of raising.
* `domain()` on a variable or parameter returns the attribute-implied
  constraints (previously an empty list).
* `matrix_frac()` evaluates an ill-conditioned second argument and warns,
  where base R's `solve()` refused and CVXPY (via NumPy) computes silently.
* `Inf` is accepted where CVXPY accepts it — constraint right-hand sides and
  variable bounds — and the decomposed `solve_via_data()` path validates
  problem data exactly like `psolve()`.
* Solver selection can skip commercial solvers:
  `options(CVXR.default_to_commercial_solvers = FALSE)` or the
  `CVXR_DEFAULT_TO_COMMERCIAL_SOLVERS` environment variable. Default
  unchanged; explicit `solver =` unaffected.

## Solver behavior

* SCS: acceleration defaults now match SCS's own C defaults, which is what
  CVXPY gets (`acceleration_lookback = 10`, interval 10); the R `scs`
  package's defaults had disabled acceleration. Numbers may move at the ~1e-5
  level; convergence on hard problems improves.
* CLARABEL: `reduced_tol_infeas_abs` is set to `5e-12` on `clarabel` 0.11.2
  and earlier (a shim — 0.11.3 fixed the default), restoring the
  almost-infeasible certificates that quasiconvex bisection needs. Pass the
  option explicitly to override.
* `accept_unknown` is now supported for CLARABEL and MOSEK, matching CVXPY:
  a stalled solve that still produced primal and dual iterates is reported as
  `optimal_inaccurate`. Deliberately opt-in; other statuses are unaffected.
* HiGHS: an infeasible problem returns a Farkas dual ray through
  `dual_value()`, as CLARABEL/OSQP/MOSEK already did. When HiGHS produces no
  ray (`ipm`/`pdlp`, bound-infeasibility), duals stay `NULL` rather than a
  zero vector that would read as a certificate. An unmapped HiGHS status is
  reported as a solver error.
* The cone scan that picks a solver looks inside constraints, not just the
  objective; and `square(x)` builds `power(x, 2)` as in CVXPY, so a squared
  term no longer reads as a power-cone requirement excluding mixed-integer
  SOCP solvers.
* A standard parameter (`feastol`/`reltol`/`abstol`/`num_iter`) with no
  mapping for the chosen solver now warns instead of being silently ignored.
* MOSEK solver parameters actually reach MOSEK now. They were merged into
  Rmosek's interface options instead of `problem$dparam`/`$iparam`, so no
  parameter ever arrived and passing one crashed with an unrelated R error.
  All three carriers work: Rmosek-native `dparam`/`iparam`/`sparam` lists,
  `MSK_*`-prefixed names, and CVXPY's `mosek_params`; unknown options raise,
  as in CVXPY. `feastol`, `reltol`, and `num_iter` translate to MOSEK's
  native tolerances; HiGHS gains `feastol` (`primal_feasibility_tolerance`).
* The `eps` keyword works as in CVXPY: for SCS it sets `eps_abs` and
  `eps_rel` (it was silently ignored); for MOSEK it fans out to all 17
  tolerance parameters, explicitly-set tolerances taking precedence.

## Performance

* About a third faster than 1.9.1 on a benchmark suite, every problem faster
  and allocating less — from 8% on the smallest problems to ~60% at 4000
  constraints. Dual recovery is now linear in the number of constraints (was
  quadratic), cutting memory up to ~60% on constraint-heavy problems.
* `geo_mean()`, `power()`, and `p_norm()` with fractional exponents are
  15-37% faster end to end.

## Notes

* PSD constraints are packed to triangular form in one `PSD -> SvecPSD`
  reduction following CVXPY 1.9, replacing per-solver copies of that logic
  (the source of the transposed-dual bug above). No user-facing API change.

# CVXR 1.9.1

This is the first CRAN release since 1.8.2 and is a large one: it tracks
CVXPY 1.9.1 and folds in the changes from the internal 1.8.2-1 and 1.9.0
development cycles. The headline additions are a derivative API for
differentiable convex programs, disciplined nonlinear programming (DNLP),
and interval-bounds propagation with native solver-bound support.

## Disciplined nonlinear programming (DNLP)

* Solve smooth nonlinear programs with `psolve(prob, nlp = TRUE)`. The
  problem must satisfy the new `is_dnlp()` grammar (disciplined
  nonlinear program); DCP problems are a subset, so any DCP problem is
  also a valid DNLP.
* New smooth atoms usable in DNLPs: `sin()`, `cos()`, `tan()`,
  `sinh()`, `tanh()`, `asinh()`, `atanh()`, and `normcdf()`.
* `prod()` is now recognized as a smooth atom, so problems involving
  products of entries (including `prod_entries(X, axis = ...)`) are
  DNLP and solvable through the NLP path.
* The classically-differentiable atoms now report `is_atom_smooth()`
  matching CVXPY 1.9: affine atoms, `exp`, `log`, `entr`, `logistic`,
  `kl_div`, `rel_entr`, `xexp`, `power` (with a constant exponent),
  `geo_mean`, `log_sum_exp`, `quad_form`, `quad_over_lin`, and `prod`.
  As a result smooth convex/concave functions such as `sum_squares()`
  are linearizable in both directions (`is_smooth()` is `TRUE`), while
  piecewise atoms (`abs`, `min`, `max`) remain convex/concave-only.
  New predicate exports: `is_smooth()`, `is_atom_smooth()`,
  `is_linearizable_convex()`, `is_linearizable_concave()`.
* Piecewise-linear atoms used in a DNLP (`abs`, `max`, `max_elemwise`,
  `norm_inf`, `sum_largest`, and `sum_smallest`) now initialize their
  epigraph variables during canonicalization, so the NLP backend has a
  complete starting point. `sum_largest` initializes its threshold to the
  (k+1)-th largest entry, matching CVXPY's warm-start fix.
* New `convolve()` atom, CVXPY 1.9's preferred (numpy-style) name for the
  1-D discrete-convolution `conv()` atom. For numeric input it falls through
  to `stats::convolve()`.
* The DNLP path is powered by derivatives from the optional
  `sparsediff` package and can solve through the optional `ipopt` and
  `Uno` R packages. NLP solver names: `"IPOPT"`, `"UNO"` (and the
  variants `"uno_ipm"` / `"uno_sqp"`), `"KNITRO"`, and `"COPT"`;
  `"KNITRO"` and `"COPT"` are registered but require solver bindings
  not yet available in R. When installed, `IPOPT` is the first automatic
  NLP solver, matching CVXPY's preference order.
* `quad_over_lin()` and `sum_squares()` now canonicalize axis-aware
  reductions (`axis = 1` and `axis = 2`) into batched second-order cone
  constraints, closing the CVXPY parity gap for
  `sum_squares(..., axis = ...)` in objectives and constraints.
* **The `UNO` interface defaults to the interior-point `ipopt` preset
  (plus MUMPS), which differs from CVXPY's `filtersqp` default.**
  CVXPY's `filtersqp` uses the BQPD active-set solver for its quadratic
  subproblems, which handles indefinite Hessians. The bundled Uno build
  is HiGHS-only (BQPD is not bundled and is not CRAN-license
  compatible), and HiGHS solves only *convex* quadratic subproblems, so
  `filtersqp` fails on nonconvex DNLPs. The interior-point `ipopt`
  preset factors the regularized KKT system with MUMPS and is robust on
  both convex and nonconvex problems. Force the SQP path with
  `solver = "uno_sqp"` (or `preset = "filtersqp"`); it is reliable only
  for problems whose subproblem Hessians stay positive semidefinite.
* Best-of-N random restarts for nonconvex DNLPs:
  `psolve(prob, nlp = TRUE, best_of = n)` solves from `n` random initial
  points and keeps the best result. Random initialization draws from a
  variable's `sample_bounds(var) <- c(low, high)` (when set) or its
  finite variable bounds. The per-run objectives are available via
  `solver_stats(prob)@extra_stats$all_objs_from_best_of`.
* After an NLP solve, `dual_value()` returns the constraint duals
  recovered from the solver (a CVXR addition; the duals are not exposed
  by CVXPY's NLP interface).

## Derivative API for differentiable convex programs

* New derivative API backed by the `diffcp` R package:
  `psolve(prob, requires_grad = TRUE)` followed by `backward(prob)` or
  `derivative(prob)`; `gradient(x)<-` and `delta(x)<-` set tangent /
  cotangent values on leaves. The chain rule is wired through `Dgp2Dcp`
  (log/exp) and `Complex2Real` (real/imag split).

## Bounds propagation

* `get_bounds()` now works on any expression, not just variables: it
  propagates interval bounds through affine, elementwise, and
  piecewise-linear atoms (e.g. `get_bounds(A %*% x + b)`,
  `get_bounds(abs(x))`, `get_bounds(sum(x))`). Bounds are returned shaped
  to the expression's dimensions.
* Variable bounds may now be sparse `Matrix` objects or symbolic bounds
  involving `Parameter`s. `get_bounds()` preserves sparse numeric bounds and
  skips symbolic bounds, while solve-time attribute lowering enforces the
  symbolic constraints. `Parameter` values containing matching `Inf` entries
  now validate correctly.
* Solvers that expose native variable bounds now consume dense numeric
  bounds directly instead of adding bound constraints: HiGHS (LP/MILP and
  QP paths), Gurobi (QP/MIQP and conic/SOCP), CPLEX (QP/MIQP), XPRESS
  (QP/MIQP and conic/SOCP), PIQP (QP), and SCIP (conic). Sparse bounds
  continue through constraint lowering.
* Parametric variable bounds use bound tensors (HiGHS), so changing
  `Parameter` bounds updates native solver bounds on DPP re-solves.
* Bounds inferred for DCP auxiliary variables are now preserved for solvers
  that consume native variable bounds.
* Attribute reduction now compacts 2D symmetric, PSD, NSD, and diagonal
  `Parameter`s before rebuilding their full expressions.
* Variables now report `Parameter`s embedded in expression bounds and include
  those bounds in DPP/DGP compliance checks.

## HiGHS warm-start and column-name validation

* HiGHS warm-start is now enabled. `psolve(prob, solver = "HIGHS",
  warm_start = TRUE)` reuses the previous solution via the persistent-solver
  API (`hi_solver_set_solution`) across LP, MILP, and QP paths. This requires
  `highs (>= 1.14)`, now available on CRAN; on older `highs` the warm-start
  tests skip and solves fall back to cold starts.
* `validate_column_name()` checks a name against the HiGHS LP-file column-name
  rules, mirroring CVXPY's `highs_conif.validate_column_name`.
* Known limitation: writing a model to a file with the original variable names
  is not yet supported — the R `highs` package exposes no column-name setter,
  so a written model would carry generic names (`c0`, `c1`, ...).

## Geometric and parameterized programming

* Positive (DGP) variables now accept numeric *and* parametric bounds under
  `gp = TRUE` (e.g. `Variable(pos = TRUE, bounds = list(lb, ub))` with `lb`/`ub`
  `Parameter`s). The DGP reduction log-transforms the bounds into the log domain
  and lowers them to constraints; parametric bounds canonicalize through the DGP
  tree without eagerly evaluating `log(value(param))`, so DPP re-solves with
  changed bound parameters work.
* `is_dpp()` gains a `context` argument (`"dcp"` or `"dgp"`), matching CVXPY's
  `is_dpp(context=...)`. In the `"dgp"` context a variable's symbolic bounds must
  be log-log-affine (e.g. a product of positive `Parameter`s is DGP-DPP though
  not DCP-DPP; `norm()`/`sum()` bounds are rejected).
* New `partial_optimize()` transform.
* CPLEX now solves LP/SOCP/MI-LP/MI-SOCP problems through the conic path by
  translating SOC blocks to `Rcplex` quadratic constraints.

## New atoms, transforms, and accessors

* New `sign()` atom (DQCP, scalar input).
* `psolve()` gains a `solver_path` argument for solver fallback chains.
* New `set_label()` / `label<-` / `format_labeled()` for attaching
  user-supplied labels to expressions.
* `power()`, `geo_mean()`, and `p_norm()` with `approx = TRUE` now
  warn when the selected solver supports power cones (`approx = FALSE`
  uses the exact `PowCone3D` / `PowConeND` form).
* MOSEK now supports mixed-integer programs (requires `Rmosek` 11.1.1+).
* `solver_stats()` for MOSEK now populates `solve_time`, `num_iters`,
  and `extra_stats`.
* New `param_dict()`, `var_dict()`, and `size_metrics()` accessors on
  `Problem`.
* `Variable(value = ...)` accepts an initial value at construction.
* New `CallbackParam` class — a `Parameter` subclass whose value is
  computed on every read by a user-supplied callback.
* `project()` accepts index-based `boolean = c(i, j, ...)` and
  `integer = c(i, j, ...)` attributes (1-based).
* `validate_arguments()` now called on `Cummax`, `Cumprod`, and
  `MinEntries` constructors.

## Performance

* Problem canonicalization and solving are faster than the 1.8.2 CRAN
  release on a wide range of problems. Two changes account for most of
  the gain: a new `.fast_new` helper for S7 expression construction
  (roughly 40% lower wall-clock on atom-dense and PSD-heavy problems),
  and routing S7 class-membership checks on the hot canonicalization
  path through a cached check on the object's class vector instead of
  `S7::S7_inherits()` (which rebuilt the class name on every call; the
  per-node cone-detection scan benefits most). Deterministic memory
  allocation is unchanged. These are internal changes with no
  user-visible API difference.

## Bug fixes

* Bug fix (also affected 1.8.x): `problem_data()` and `get_problem_data()` now
  take an explicit `gp` argument. Previously `gp = TRUE` passed through `...` was
  silently ignored, compiling a geometric program as a DCP problem.
* Bug fix (also affected 1.8.x): `psolve()` now takes explicit `enforce_dpp` and
  `ignore_dpp` arguments, matching CVXPY's `solve()`. Previously they were
  silently swallowed by `...`; `enforce_dpp = TRUE` now raises on a non-DPP
  parametrized problem instead of falling back to a non-DPP compile.
* Fixed QP-path canonicalization of `quad_over_lin(expr, constant)` so the
  decision uses the canonicalized denominator, matching CVXPY 1.9 behavior for
  cases such as `quad_over_lin(exp(x), 1)` that should keep `ExpCone` without
  adding an SOC block.
* Fix `Variable(c(n, n), diag = TRUE)` handling through `CvxAttr2Constr`.
* Fix `perspective()` canonicalizer on PSD, NSD, and diag matrix variables.
* Fix `project()` on sparse `Matrix`-package objects.

## Notes

* Complex-parameter DPP fast-path parity remains a deliberate feature gap:
  complex parameters are evaluated before `Complex2Real` because CVXR's R
  sparse-matrix backend cannot carry complex sparse parameter tensors. Ordinary
  complex expression solving is still supported through `Complex2Real`.
* CVXR's 2-D positive-dimensional model does not support zero-sized
  expressions, constraints, or solves. Generated or data-dependent models
  should skip vacuous constraints and use explicit zero contributions when a
  mask or filter selects no entries.

# CVXR 1.8.2

## New solvers: SCIP and XPRESS (15 total)

* Added SCIP solver (14th solver) via the R `scip` package. Supports
  LP, SOCP, MI-LP, and MI-SOCP problems. SCIP is registered in the
  conic solver path with `SUPPORTED_CONSTRAINTS = list(Zero, NonNeg,
  SOC)` and `MIP_CAPABLE = TRUE`.
* SCIP solver parameters can be passed via `scip_params` sub-list
  (matching CVXPY convention) for path-style parameter names like
  `"limits/time"`, `"limits/gap"`.
* Duals are not currently extracted (R `scip` package lacks dual API).
* 29 tests mirroring CVXPY's `TestSCIP` class.
* Added XPRESS solver (15th solver) via the R `xpress` package with
  dual-interface architecture: QP path (`XPRESS_QP_Solver`) for LP/QP
  and conic path (`XPRESS_Conic_Solver`) for SOCP/MI-LP/MI-SOCP.
  Both paths support Zero and NonNeg constraints; the conic path also
  supports SOC. MIP warm-start is supported.

## `diag()` and `norm()` dispatch fixes

* `diag()` now works on CVXR expressions: `diag(vector_expr)` creates
  a diagonal matrix (`DiagVec`), `diag(square_matrix_expr)` extracts
  the diagonal (`DiagMat`). Falls through to `Matrix::diag()` for
  non-Expression inputs, correctly handling both Matrix S4 objects
  and base R matrices.
* `norm()` fallthrough now delegates to `Matrix::norm()` instead of
  `base::norm()`, fixing `norm(sparse_matrix)` when CVXR is loaded.
* `t()` and `mean()` switched from S7 `method()` to
  `registerS3method()` to avoid demoting Matrix's S4 generics.

## CVXPY 1.8.2 parity

All applicable bug fixes from CVXPY v1.8.2 have been ported, each
annotated with `## CVXPY v1.8.2 fix:` in the source.

## `solver_opts()` and `use_quad_obj`

* New `solver_opts()` constructor for unified solver options. Standard
  tolerance parameters (`feastol`, `reltol`, `abstol`, `num_iter`) and
  solver-specific parameters now flow through `psolve(...)` via the
  `solver_opts()` constructor.
* New `use_quad_obj` option (default `TRUE`): when `FALSE`, forces
  conic decomposition path instead of QP, enabling `quad_form_canon`
  to detect indefinite P matrices.
* New `supports_quad_obj()` generic on conic solvers (Clarabel and SCS
  return `TRUE`, others `FALSE`).

## Element-wise matrix indexing (`SpecialIndex`)

* New `SpecialIndex` atom (mirroring CVXPY's `special_index`) enables
  R-idiomatic element-wise selection on matrix expressions:
    - 2-column matrix indexing: `x[cbind(rows, cols)]`
    - Logical matrix indexing: `x[mask]`
    - Linear integer indexing: `x[c(1, 5, 9)]`
    - Logical vector indexing: `x[c(TRUE, FALSE, ...)]`
* This makes it natural to constrain specific entries of a matrix
  variable, e.g., fixing observed entries of a partially-known matrix:
  ```r
  ind <- which(!is.na(Rmiss), arr.ind = TRUE)
  prob <- Problem(Minimize(obj), list(X[ind] == Rmiss[ind]))
  ```
* Previously, `x[i]` on a matrix expression errored with "Single-index
  selection on matrices not supported." The workaround required an
  element-wise multiply with a binary mask, which was unintuitive.
* Internally uses sparse selection matrix multiplication (reshape →
  sparse matmul), requiring no C++ changes.

## Bug Fixes

* Fixed O(n²) memory usage in `sum_squares()`, `power(x, 2)`,
  `quad_over_lin()`, and `huber()` when the first argument has many
  elements (e.g., regression residuals with thousands of observations).
  The quadratic canonicalizer was converting a sparse identity matrix
  to dense form. Memory now scales linearly with problem size.

# CVXR 1.8.1

## Complete rewrite using S7 object system

This is a ground-up rewrite of CVXR using R's S7 object system,
designed to be isomorphic with CVXPY 1.8.1 for long-term
maintainability. ~4-5x faster than CVXR 1.0-15 on typical problems.

### New features

* S7 class system replaces S4 for all expression, constraint, and
  problem classes. Significantly faster construction and method
  dispatch.
* 13 solvers at initial release: Clarabel (default), SCS, OSQP, HiGHS,
  MOSEK, Gurobi, GLPK, GLPK_MI, ECOS, ECOS_BB, CPLEX, CVXOPT, PIQP.
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
* Full 2D broadcasting for `+` and `-` operators, matching CVXPY's
  `Expression.broadcast()`. Expressions with compatible but different
  shapes (e.g., `(1, n) + (m, 1)` → `(m, n)`) now work correctly in
  constraints and objectives. Previously only scalar promotion was
  supported.
* Correct solver-cone routing matching CVXPY's atom classifications.
  Atoms that canonicalize to PSD cones (`MatrixFrac`, `TrInv`,
  `SigmaMax`, `NormNuc`, `LambdaSumLargest`) are now correctly
  classified, ensuring solvers like ECOS that don't support SDP are
  rejected with a clear error message instead of a cryptic dimension
  mismatch. Approximate variants (`PnormApprox`, `PowerApprox`,
  `GeoMeanApprox`) correctly route to SOC; exact variants route to
  their native cone types (`PowCone3D`, `PowConeND`).

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
* Warm-start for HiGHS is implemented (works with `highs` 1.14) but the
  upstream R `highs` PR exposing `setSolution()` has not yet been
  actioned by the maintainers; lives on branch `highs-warm-start`.
* Complex SOC dual variables are not recovered (CVXPY itself does not
  implement this — `complex2real.py:56` declares
  `UNIMPLEMENTED_COMPLEX_DUALS = (SOC, OpRelEntrConeQuad)`). Complex
  SOC *primal* canonicalization works.
* Deferred to a future release: RelEntrConeQuad / OpRelEntrConeQuad,
  quantum atoms, CPLEX conic (SOC) path.

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
