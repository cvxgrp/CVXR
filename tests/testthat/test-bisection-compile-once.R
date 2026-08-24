## The compile-once fast path in .lower_problem (ADR D_PERF.9).
##
## CVXPY SOURCE: bisection.py:29-31 reconstructs a fresh Problem for EVERY
## bisection query, and bisection.py:47-52 is upstream's own TODO wishing it
## did not: each reconstruction forces a full canonicalization, so a DQCP
## solve pays ~20-50 compiles for one answer.  When Dqcp2Dcp emits no lazy
## constraints, the parameterized problem IS its own lowering and is DPP
## (the bisection parameter enters as param * variable), so CVXR returns the
## SAME object and repeated solves ride the parameterized-problem cache.
## Measured on Maximize(sqrt(x)/exp(x)) with SCS: 1.45s -> 0.10s.
##
## This is a deliberate, recorded deviation from upstream -- see
## notes/decisions.md D_PERF.9.

## @cvxpy NONE -- R-specific deviation D_PERF.9 (upstream TODO, unimplemented
## there); guards that the fast path actually engages.
test_that(".lower_problem returns the problem itself when nothing is lazy", {
  x <- Variable()
  problem <- Problem(Maximize(sqrt(x) / exp(x)))
  chain <- Chain(reductions = list(FlipObjective(), Dqcp2Dcp()))
  param_problem <- reduction_apply(chain, problem)[[1L]]

  expect_length(param_problem@.cache$lazy_constraints, 0L)
  expect_true(is_dpp(param_problem))
  ## Pre-D_PERF.9 this was a fresh Problem every call; identity is what lets
  ## the parameterized-problem cache survive across bisection iterations.
  expect_identical(.lower_problem(param_problem), param_problem)
})

## @cvxpy NONE -- the load-bearing branch: with lazy (sign-branching)
## constraints the per-query reconstruction must be kept.
test_that(".lower_problem still reconstructs when lazy constraints exist", {
  x <- Variable(5)
  problem <- Problem(Minimize(length_expr(x)), list(sum(x) >= 1))
  param_problem <- reduction_apply(Chain(reductions = list(Dqcp2Dcp())),
                                   problem)[[1L]]

  expect_length(param_problem@.cache$lazy_constraints, 1L)
  t_param <- param_problem@.cache$bisection_data$param
  value(t_param) <- 5
  lowered <- .lower_problem(param_problem)
  expect_false(identical(lowered, param_problem))
  ## The lazy callable materialized: one more constraint than the base set.
  expect_length(lowered@constraints,
                length(param_problem@constraints) + 1L)
})

## @cvxpy test_dqcp.py::TestDqcp::test_concave_frac -- end-to-end through the
## fast path; the docs example that motivated D_PERF.9.
test_that("concave-fractional solves through the fast path", {
  x <- Variable()
  problem <- Problem(Maximize(sqrt(x) / exp(x)))
  expect_true(is_dqcp(problem))
  ## SCS, not the default: near the optimal ratio the feasibility
  ## subproblems are marginally infeasible, and CLARABEL/MOSEK end in
  ## NumericalError there instead of certifying (CVXPY 1.9.2 fails
  ## identically with its default).  SCS and ECOS certify cleanly.
  val <- suppressWarnings(psolve(problem, qcp = TRUE, solver = "SCS"))
  expect_true(status(problem) %in% SOLUTION_PRESENT)
  expect_equal(val, sqrt(0.5) * exp(-0.5), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy NONE -- fallback end-to-end: a lazy-constraint DQCP problem must
## still solve exactly as before the fast path existed.
test_that("lazy-constraint DQCP problems still solve", {
  x <- Variable(5)
  problem <- Problem(Minimize(length_expr(x)), list(sum(x) >= 1))
  val <- suppressWarnings(psolve(problem, qcp = TRUE, solver = "SCS"))
  expect_equal(val, 1)
  expect_true(status(problem) %in% SOLUTION_PRESENT)
})
