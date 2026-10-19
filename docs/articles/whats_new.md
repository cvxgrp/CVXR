# What's New in CVXR

This vignette highlights notable user-facing changes in each release of
the S7 rewrite of `CVXR`, newest first. For the complete, fine-grained
list see the package `NEWS` file (`news(package = "CVXR")`). For the
introductory tutorial, see
[`vignette("cvxr_intro")`](https://www.cvxgrp.org/CVXR/articles/cvxr_intro.md);
for worked examples, visit the [CVXR website](https://cvxr.rbind.io).

- [CVXR 1.9.3](#cvxr-193) — correctness fixes, matrix norms as in base
  R, CVXPY 1.9 solver selection
- [CVXR 1.9.2](#cvxr-192) — correctness fixes, solver-options overhaul,
  faster solves
- [CVXR 1.9.1](#cvxr-19x) — disciplined nonlinear programming,
  derivatives, bounds propagation, new atoms
- [CVXR 1.8.x](#cvxr-18x) — the ground-up S7 rewrite

## CVXR 1.9.3

`CVXR` 1.9.3 tracks CVXPY 1.9.3. Like 1.9.2, it is chiefly a
**correctness release**: several defects that returned wrong answers
with no error or warning are fixed. Matrix norms now mean what they mean
in base R, and the solver is chosen the way CVXPY 1.9 chooses it.

### Correctness

These fixes change results. The new results are the answers to the
problem you wrote.

- **Gurobi** solved quadratic objectives with cross terms, such as
  `quad_form(w, Sigma)` in a portfolio problem, with every cross term
  halved. This was present in 1.9.2. Gurobi’s second-order-cone duals
  also had the wrong sign.
- **CVXOPT** reported feasible problems with only equality constraints
  as infeasible, and unbounded problems as optimal.
- [`perspective()`](https://www.cvxgrp.org/CVXR/reference/perspective.md)
  of a function of a `symmetric`, `PSD` or `NSD` matrix variable
  returned a wrong optimum.
- A DGP problem with parameters could be reported DPP and solved with
  wrong solver data.
- [`gmatmul()`](https://www.cvxgrp.org/CVXR/reference/gmatmul.md) with
  parameters on both sides gave a wrong answer under `gp = TRUE`.
- `psolve(solver = "CLARABEL", warm_start = TRUE)` returned the previous
  optimum when the cones had changed.
- A complex variable declared `PSD` or `NSD` lost its cone.
- `sum_entries(x, axis = 2)` of a one-row expression gave wrong answers
  or crashed R. Present in 1.9.2.
- `%*%` with a scalar operand that does not conform, as in `2 %*% x` for
  a 2 x 3 `x`, or with a 1 x 1 operand that does not conform, built a
  product that gave a scalar product or crashed R. Both are now errors,
  as in base R; use `*` to multiply by a scalar. A scalar operand that
  conforms, as in `x %*% 2` for a 3 x 1 `x`, is a product, as in base R.
- [`sum()`](https://rdrr.io/r/base/sum.html),
  [`max()`](https://rdrr.io/r/base/Extremes.html),
  [`min()`](https://rdrr.io/r/base/Extremes.html) and
  [`prod()`](https://rdrr.io/r/base/prod.html) of an expression took a
  named `axis` as one more operand, so `sum(x, axis = 2)` was
  `sum(x + 2)`. They now honor `axis` and `keepdims`.
- [`sum()`](https://rdrr.io/r/base/sum.html),
  [`prod()`](https://rdrr.io/r/base/prod.html),
  [`max()`](https://rdrr.io/r/base/Extremes.html) and
  [`min()`](https://rdrr.io/r/base/Extremes.html) with more than one
  operand now mean what they mean in base R: `sum(x, 2)` is `sum(x) + 2`
  and `max(x, y)` is one scalar, the largest entry of either. They were
  `sum(x + 2)`, `prod(x * y)` and the elementwise maximum and minimum;
  use
  [`max_elemwise()`](https://www.cvxgrp.org/CVXR/reference/max_elemwise.md)
  and
  [`min_elemwise()`](https://www.cvxgrp.org/CVXR/reference/min_elemwise.md)
  for those.
- [`max_elemwise()`](https://www.cvxgrp.org/CVXR/reference/max_elemwise.md)
  and
  [`min_elemwise()`](https://www.cvxgrp.org/CVXR/reference/min_elemwise.md)
  take `na.rm`, and
  [`hstack()`](https://www.cvxgrp.org/CVXR/reference/hstack.md) and
  [`vstack()`](https://www.cvxgrp.org/CVXR/reference/vstack.md) take
  `deparse.level`, as [`pmax()`](https://rdrr.io/r/base/Extremes.html)
  and [`rbind()`](https://rdrr.io/r/base/cbind.html) do, with no effect.
  `max_elemwise(u, v, na.rm = TRUE)` used to take `TRUE` as a third
  operand. Named operands, as from `do.call(vstack, split(...))`, are
  accepted as before.
- `x / A` with a matrix or vector `A` divided every entry of `x` by
  `A[1, 1]`. It is now entrywise, as in CVXPY. Present in 1.9.2.
- `x / d` with shapes that do not broadcast, such as a 2 x 3 `x` and a
  length-3 `d`, recycled the divisor silently. It is now an error, as in
  CVXPY; a row or column divisor still broadcasts.
- [`value()`](https://www.cvxgrp.org/CVXR/reference/value.md) of a sum,
  a product or an index of a sparse-matrix constant failed with
  “non-conformable arrays” or “subscript out of bounds”, and a sum of
  two sparse constants came back flattened to one column. These now
  evaluate; a mixed sparse and dense result is a base matrix, a sparse
  one stays sparse.
- Constants from Matrix classes other than `dgCMatrix` were misread: a
  unit `Diagonal(n)` was all zeros, so `quad_form(x, Diagonal(n))` was
  unbounded. Present in 1.9.2.
- A sparse `bounds` matrix bound only its stored entries. As in CVXPY,
  it is now accepted only on a `sparsity` variable with the same
  pattern. Present in 1.9.2.
- A CLARABEL or OSQP warm re-solve after a parameter became zero, or
  nonzero, returned the previous optimum as optimal. Present in 1.9.2.
- [`conv()`](https://www.cvxgrp.org/CVXR/reference/conv.md),
  [`cumsum()`](https://rdrr.io/r/base/cumsum.html) and
  [`kron()`](https://www.cvxgrp.org/CVXR/reference/kron.md) were treated
  as log-log affine, so a geometric program containing them was accepted
  and could solve to a wrong value. Such a problem is now rejected, as
  in CVXPY.
- Calling
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
  on a geometric program between two solves made the next solve after a
  parameter change return the previous optimum.
  [`problem_data()`](https://www.cvxgrp.org/CVXR/reference/problem_data.md)
  now goes through the same compile cache as
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md).
- [`max_elemwise()`](https://www.cvxgrp.org/CVXR/reference/max_elemwise.md)
  and
  [`min_elemwise()`](https://www.cvxgrp.org/CVXR/reference/min_elemwise.md)
  of a column and a row recycled instead of broadcasting, so
  [`psolve()`](https://www.cvxgrp.org/CVXR/reference/psolve.md) returned
  a wrong objective.
- On an expression, `log(x, base)` ignored `base`, `outer(x, y, FUN)`
  ignored `FUN` and `diag(x, nrow)` ignored `nrow`. `base` is now
  honored; a `FUN` other than `"*"`, or a `nrow`, is an error.
- `x[NULL, ]` of an expression returned all of `x`, and
  `x[c(TRUE, NA, FALSE), ]` dropped the `NA`. Both are now errors.
- `cumsum(X)` of a matrix expression is now column-wise, as in CVXPY and
  as [`cumprod()`](https://rdrr.io/r/base/cumsum.html) and
  [`cummax()`](https://rdrr.io/r/base/cumsum.html) already were.

### Matrix norms follow base R

With no `axis`, `cvxr_norm(X, p)` and `norm(X, type)` on a matrix now
agree with base R’s
[`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) and with
CVXPY:

- `p = 1` is the maximum absolute column sum.
- `p = 2` is the largest singular value.
- `p = Inf` is the maximum absolute row sum.

They used to be entrywise. A row vector is a 1 x n matrix, as in base R.
Any other `p` on a matrix is now an error. For the entrywise norms use
`cvxr_norm(X, "fro")`, `norm1(X)` or `p_norm(X, p)`.
`cvxr_norm(X, "nuc")` now gives the nuclear norm. With no `type`,
`norm(X)` of an expression is the spectral norm, as in CVXPY; on an
ordinary R matrix it is base R’s one-norm.

\
`X`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``-``2``, ``3``, ``4``)``, ``2``, ``2``)`\
[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(`[`cvxr_norm`](https://www.cvxgrp.org/CVXR/reference/cvxr_norm.md)`(`[`Constant`](https://www.cvxgrp.org/CVXR/reference/Constant.md)`(``X``)``, ``1``)``)``   ``# 7, as norm(X, "1")`\
[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(`[`norm1`](https://www.cvxgrp.org/CVXR/reference/norm1.md)`(`[`Constant`](https://www.cvxgrp.org/CVXR/reference/Constant.md)`(``X``)``)``)``          ``# 10, the entrywise 1-norm`

### Base R functions keep their meaning

`CVXR` masks
[`diag()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md),
[`norm()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md),
[`outer()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md),
[`sd()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md),
[`var()`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md) and
[`convolve()`](https://www.cvxgrp.org/CVXR/reference/convolve.md) so
that they also work on expressions. On anything that is not an
expression, each is now exactly the base R or stats function, called as
written. On an expression each follows CVXPY, and an argument meant only
for ordinary data, such as `na.rm`, is an error.

[`power()`](https://www.cvxgrp.org/CVXR/reference/cvxr_power.md) is
different: it is CVXPY’s `power(x, p)` for any argument, as before, and
the new
[`cvxr_power()`](https://www.cvxgrp.org/CVXR/reference/cvxr_power.md) is
the same function. It masks the unrelated
[`stats::power()`](https://rdrr.io/r/stats/power.html), a
[`glm()`](https://rdrr.io/r/stats/glm.html) link, so with `CVXR`
attached write [`stats::power()`](https://rdrr.io/r/stats/power.html) in
a [`glm()`](https://rdrr.io/r/stats/glm.html) family.

### Solver selection

The solver is now chosen as in CVXPY 1.9: first from the problem’s
structure, and then the solver’s QP or conic interface follows. As a
result, a problem may be solved by a different installed solver than in
1.9.2:

- Commercial solvers are preferred for QPs and mixed-integer QPs too.
- SCIP solves mixed-integer QPs.
- A named HiGHS, Gurobi, CPLEX or Xpress solves an LP through its conic
  interface.
- Falling back to a non-default solver warns.
- ECOS_BB is never chosen unless you name it.

When no installed solver fits, the error reads “No installed solver
could handle this problem.”

### Solvers

- **MOSEK** takes its status from the solution status, as CVXPY does, so
  an iteration or time limit is no longer an error by itself. The
  simplex optimizers work.
- **HiGHS** accepts options nested in a `highs_options` list. An option
  HiGHS rejects is now an error rather than being dropped.
- **Gurobi** `warm_start = TRUE` seeds the start point from the
  variables’ values.
- **CPLEX** returns dual values for continuous problems with
  second-order cones; they were all missing.
- **SCS** exposes the infeasibility certificate of an infeasible solve
  through
  [`dual_value()`](https://www.cvxgrp.org/CVXR/reference/dual_value.md).
- `solver_stats(prob)@extra_stats` now carries the raw solver result for
  OSQP, ECOS, Gurobi, HiGHS and PIQP as well.
- **GLPK_MI** reports an infeasible mixed-integer problem as
  `infeasible`; it was a solver error.
- Solving no longer advances R’s random number stream, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) followed by a solve
  and then random draws gives the same draws as without the solve.
  Present in every earlier version.

### New features

- [`geo_mean()`](https://www.cvxgrp.org/CVXR/reference/geo_mean.md)
  accepts arguments of any shape and gains `axis` and `keepdims`.
  [`cvxr_var()`](https://www.cvxgrp.org/CVXR/reference/cvxr_var.md)
  gains them too.
- Exact `geo_mean(x, approx = FALSE)` and `PowConeND` constraints now
  solve on SCS and MOSEK as well as Clarabel.
- [`curvature()`](https://www.cvxgrp.org/CVXR/reference/curvature.md)
  reports log-log and quasi curvature instead of `"UNKNOWN"`. The new
  [`curvatures()`](https://www.cvxgrp.org/CVXR/reference/curvatures.md)
  returns every curvature that applies, and
  [`log_log_curvature()`](https://www.cvxgrp.org/CVXR/reference/log_log_curvature.md)
  and [`is_diag()`](https://www.cvxgrp.org/CVXR/reference/is_diag.md)
  are new as well.
- A quadratic objective over a structured variable (`symmetric`, `PSD`,
  `diag`, …) now solves on QP solvers.

\
`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``pos ``=`` ``TRUE``)`\
`y`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``pos ``=`` ``TRUE``)`\
[`curvature`](https://www.cvxgrp.org/CVXR/reference/curvature.md)`(``x`` ``*`` ``y``)``                ``# "LOG-LOG AFFINE"`\
[`curvature`](https://www.cvxgrp.org/CVXR/reference/curvature.md)`(`[`ceiling`](https://rdrr.io/r/base/Round.html)`(`[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``)``)``)``  ``# "QUASILINEAR"`

### Performance

- The first solve of a parametric (DPP) problem compiles about twice as
  fast and uses far less memory.
- A DPP problem with a parametric quadratic objective keeps its compiled
  form across solves.
- `psolve(nlp = TRUE)` no longer makes dense copies of sparse constant
  matrices.

## CVXR 1.9.2

`CVXR` 1.9.2 tracks CVXPY 1.9.2. It is primarily a **correctness
release**: several defects that returned plausible-looking wrong answers
with no error or warning are fixed, solver options are easier to pass
and harder to lose silently, and solving is faster across the board.

### Correctness: weighted quadratic objectives

A constant weight multiplying `sum(square(x))` was discarded when the
quadratic objective was assembled, so penalized least-squares problems
were solved with the wrong penalty — ridge regression being the common
case, where every value of the penalty parameter returned the same fit.
The defect had been present since 1.8.2;
[`sum_squares()`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md),
[`quad_form()`](https://www.cvxgrp.org/CVXR/reference/quad_form.md), and
`p_norm(x, 2)^2` were never affected. If you used a weighted
`sum(square(x))` objective, results can legitimately change with this
release — they are now the answers to the problem you wrote.

### `dual_value()` returns constraint-shaped duals

[`dual_value()`](https://www.cvxgrp.org/CVXR/reference/dual_value.md)
now returns a matrix shaped like its constraint rather than a flat
vector, making the documented return type accurate. Linear indexing
(`d[3]`) is unchanged, since R matrices index linearly in column-major
order, but code calling [`length()`](https://rdrr.io/r/base/length.html)
or [`dim()`](https://rdrr.io/r/base/dim.html) on a dual may need
adjusting. Scalar constraints, `SOC`, and `ExpCone` duals are
deliberately unchanged, matching CVXPY.

### Solver options

Passing solver options is now more uniform, and options can no longer
vanish silently:

- The standard tolerance names `feastol`, `reltol`, `abstol`, and
  `num_iter` now also map onto MOSEK and HiGHS.
- A standard option that a solver cannot express now produces a
  **warning** (condition class `cvxr_unmapped_solver_param`) instead of
  being silently dropped.
- **MOSEK** accepts CVXPY-style `mosek_params` and the `eps` keyword,
  which fans out to MOSEK’s tolerance parameters exactly as in CVXPY;
  all options are routed to the correct `Rmosek` parameter groups.
- **SCS** accepts `eps`, which sets both `eps_abs` and `eps_rel`.

\
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``, solver ``=`` ``"MOSEK"``, eps ``=`` ``1e-8``)`\
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``, solver ``=`` ``"SCS"``, eps ``=`` ``1e-6``)`

### New warnings and catchable errors

- A parameterized problem that is not DPP now warns that repeated solves
  will not be faster, matching CVXPY; the answers are unchanged.
- [`matrix_frac()`](https://www.cvxgrp.org/CVXR/reference/matrix_frac.md)
  with an ill-conditioned second argument now returns a value (as CVXPY
  does) and warns, where it previously raised R’s “computationally
  singular” error.
- Five error conditions — `DPPError`, `DGPError`, `DQCPError`,
  `DNLPError`, and `ParameterError` — are now catchable by class with
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html), as in CVXPY;
  previously they could only be matched by message text.

### Mixed-integer declarations

A *partial* `boolean`/`integer` index list (for example
`Variable(2, boolean = c(1))`) now correctly makes the problem
mixed-integer. Such problems previously appeared to solve with a
continuous solver; they now require a mixed-integer-capable solver and
say so if none is installed.

### Performance

1.9.2 is roughly a third faster than 1.9.1 overall, with substantially
lower memory use on larger problems, and DQCP problems solve
dramatically faster: the bisection now compiles its subproblem once and
re-solves it per iterate instead of rebuilding it from scratch.

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

\
`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``2``)`\
`prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(`[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum_squares`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)`(``x`` ``-`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2``)``)``)``)`\
[`is_dnlp`](https://www.cvxgrp.org/CVXR/reference/is_dnlp.md)`(``prob``)``                       ``# TRUE`\
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``, nlp ``=`` ``TRUE``)``            ``# solved through the NLP path`

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

\
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, requires_grad ``=`` ``TRUE``)`\
[`delta`](https://www.cvxgrp.org/CVXR/reference/delta.md)`(``a``)`` ``<-`` ``da``                      ``# perturbation of parameter a`\
[`derivative`](https://www.cvxgrp.org/CVXR/reference/derivative.md)`(``problem``)``                 ``# propagate forward`\
[`delta`](https://www.cvxgrp.org/CVXR/reference/delta.md)`(``x``)``                            ``# resulting change in variable x`

**Reverse mode** (gradient of the solution with respect to parameters):

\
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``problem``, requires_grad ``=`` ``TRUE``)`\
[`backward`](https://www.cvxgrp.org/CVXR/reference/backward.md)`(``problem``)``                   ``# propagate backward`\
[`gradient`](https://www.cvxgrp.org/CVXR/reference/gradient.md)`(``a``)``                         ``# d(solution) / d(a)`

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

\
`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``3``, bounds ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``-``1``, ``2``)``)`\
[`get_bounds`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)`(``A`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``x`` ``+`` ``b``)``         ``# bounds propagated through the affine map`\
[`get_bounds`](https://www.cvxgrp.org/CVXR/reference/get_bounds.md)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(``x``)``)``             ``# and through atoms`

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

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`CVXR`](https://cvxr.rbind.io)`)`\
`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``2``)`\
`prob`` ``<-`` `[`Problem`](https://www.cvxgrp.org/CVXR/reference/Problem.md)`(`[`Minimize`](https://www.cvxgrp.org/CVXR/reference/Minimize.md)`(`[`sum_squares`](https://www.cvxgrp.org/CVXR/reference/sum_squares.md)`(``x``)``)``, `[`list`](https://rdrr.io/r/base/list.html)`(``x`` ``>=`` ``1``)``)`\
`opt_val`` ``<-`` `[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)`(``prob``)``       ``# returns optimal value directly`\
`x_val`` ``<-`` `[`value`](https://www.cvxgrp.org/CVXR/reference/value.md)`(``x``)``             ``# extract variable value`\
`prob_status`` ``<-`` `[`status`](https://www.cvxgrp.org/CVXR/reference/status.md)`(``prob``)``   ``# check status`

The old [`solve()`](https://rdrr.io/r/base/solve.html) still works but
returns a backward-compatible list:

\
`result`` ``<-`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``prob``)`\
`result``$``value``       ``# optimal value`\
`result``$``getValue``(``x``)`` ``# variable value (deprecated)`\
`result``$``status``      ``# problem status`

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

The current table of solvers and the problem classes each one handles is
kept on the website: [Solver
Capabilities](https://cvxr.rbind.io/examples/solvers/using-other-solvers.html#solver-capabilities).

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

\
`x`` ``<-`` `[`Variable`](https://www.cvxgrp.org/CVXR/reference/Variable.md)`(``3``)`\
[`abs`](https://rdrr.io/r/base/MathFun.html)`(``x``)``        ``# elementwise absolute value`\
[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``x``)``       ``# elementwise square root`\
[`sum`](https://rdrr.io/r/base/sum.html)`(``x``)``        ``# sum of entries`\
[`max`](https://rdrr.io/r/base/Extremes.html)`(``x``)``        ``# maximum entry`\
[`norm`](https://www.cvxgrp.org/CVXR/reference/math_atoms.md)`(``x``, ``"2"``)``  ``# Euclidean norm`

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

    \
    `## Old (worked in CVXR 1.x because axis reductions were 1D):`\
    [`sum_entries`](https://www.cvxgrp.org/CVXR/reference/sum_entries.md)`(``X``, axis ``=`` ``2``)`` ``==`` ``target_vec`\
    `## New (wrap target as row vector to match the (1, n) shape):`\
    [`sum_entries`](https://www.cvxgrp.org/CVXR/reference/sum_entries.md)`(``X``, axis ``=`` ``2``)`` ``==`` `[`t`](https://rdrr.io/r/base/t.html)`(``target_vec``)`

    Similarly, if you extract a scalar from a CVXR result and need a
    plain numeric value, use
    [`as.numeric()`](https://rdrr.io/r/base/numeric.html) to drop the
    matrix dimensions.

### CRAN Submission Tip

If you encounter issues involving the `Rmosek` package while submitting
your package to CRAN, include the following code in `<your_pkg>/R/zzz.R`
to resolve the issue.

\
`## Content of <your_pkg>/R/zzz.R`\
\
`.onLoad`` ``<-`` ``function``(``libname``, ``pkgname``)`` ``{`\
`  ``CVXR``::`[`exclude_solvers`](https://www.cvxgrp.org/CVXR/reference/available_solvers.md)`(``"MOSEK"``)`\
`}`\
`.onUnload`` ``<-`` ``function``(``libname``, ``pkgname``)`` ``{`\
`  ``CVXR``::`[`include_solvers`](https://www.cvxgrp.org/CVXR/reference/available_solvers.md)`(``"MOSEK"``)`\
`}`

### Further Reading

- [CVXR website](https://cvxr.rbind.io) — worked examples
- [Package reference](https://www.cvxgrp.org/CVXR/) — full API
  documentation
- [CVXPY documentation](https://www.cvxpy.org/) — mathematical framework
- Fu, Narasimhan, and Boyd (2020). “CVXR: An R Package for Disciplined
  Convex Optimization.” *Journal of Statistical Software*, 94(14).
  <doi:10.18637/jss.v094.i14>
