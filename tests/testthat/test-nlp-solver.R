# Tests for nlp_solver: NLPsolver + Bounds + Oracles (CVXPY 1.9.0, Step 5b).
#
# NLPsolver.apply builds NLP problem data; Bounds lowers a problem to the NLP
# standard form (g^l <= g(x) <= g^u, x^l <= x <= x^u); Oracles wraps the
# diff_engine C_problem for value/derivative evaluation. The runnable solve
# (Uno) is uno_nlpif, now landed and active.
#
# CVXPY SOURCE: reductions/solvers/nlp_solvers/nlp_solver.py

library(CVXR)

# ===================================================================
# Bounds: lowering + constraint/variable bounds + initial point
# ===================================================================

## @cvxpy NONE
test_that("Bounds lowers constraints and extracts cl/cu/lb/ub/x0", {
  x <- Variable(2, nonneg = TRUE)
  value(x) <- c(1, 2)
  prob <- Problem(Minimize(sum_entries(Exp(x))),
                  list(sum_entries(x) == 3, sum_entries(x) <= 10))
  b <- Bounds(prob)
  ## equality -> [0,0]; inequality -> [0, Inf]
  expect_equal(b@cl, c(0, 0))
  expect_equal(b@cu, c(0, Inf))
  ## nonneg variable bounds
  expect_equal(b@lb, c(0, 0))
  expect_equal(b@ub, c(Inf, Inf))
  ## initial point = variable value
  expect_equal(b@x0, c(1, 2))
  ## lowered constraint types: Equality -> Zero, Inequality -> NonNeg
  expect_true(S7_inherits(b@new_problem@constraints[[1L]], Zero))
  expect_true(S7_inherits(b@new_problem@constraints[[2L]], NonNeg))
})

## @cvxpy NONE
test_that("Bounds: a plain variable has (-Inf, Inf) bounds", {
  x <- Variable(3)
  value(x) <- c(0, 0.5, 1)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
  b <- Bounds(prob)
  expect_equal(b@lb, rep(-Inf, 3))
  expect_equal(b@ub, rep(Inf, 3))
})

## @cvxpy NONE
test_that("Bounds aborts when a variable has no initial value", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list())
  expect_error(Bounds(prob), "no value")
})

# ===================================================================
# NLPsolver: accepts + apply
# ===================================================================

## @cvxpy NONE
test_that("NLPsolver accepts a DNLP problem and apply yields NLP data", {
  x <- Variable(3)
  value(x) <- c(0, 0.5, 1)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
  solver <- NLPsolver()
  expect_true(solver@BOUNDED_VARIABLES)
  expect_true(reduction_accepts(solver, prob))
  res <- reduction_apply(solver, prob)
  data <- res[[1L]]
  expect_equal(data$cl, c(0))
  expect_equal(data$cu, c(0))
  expect_equal(data$x0, c(0, 0.5, 1))
  expect_true(S7_inherits(data[["_bounds"]], Bounds))
})

# ===================================================================
# Oracles: value/gradient/Jacobian/Hessian through the diff engine
# ===================================================================

## @cvxpy NONE
test_that("Oracles wrap C_problem and evaluate sum(exp(x)) derivatives", {
  skip_if_not_installed("sparsediff")
  x <- Variable(3)
  value(x) <- c(0, 0.5, 1)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
  b <- Bounds(prob)
  orc <- .nlp_oracles(b@new_problem, use_hessian = TRUE)
  u <- c(0, 0.5, 1)
  expect_equal(orc$objective(u), sum(exp(u)), tolerance = 1e-6)
  expect_equal(as.numeric(orc$gradient(u)), exp(u), tolerance = 1e-6)
  orc$constraints(u)
  expect_equal(as.numeric(orc$jacobian(u)), c(1, 1, 1), tolerance = 1e-9)
  h <- as.numeric(orc$hessian(u, duals = c(0), obj_factor = 1))
  expect_equal(sort(h), sort(exp(u)), tolerance = 1e-6)
  ## sparsity structures are available and cached
  expect_false(is.null(orc$jacobianstructure()))
  expect_false(is.null(orc$hessianstructure()))
})

## @cvxpy NONE
test_that("Oracles with use_hessian = FALSE returns empty Hessian structure", {
  skip_if_not_installed("sparsediff")
  x <- Variable(2)
  value(x) <- c(1, 1)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list())
  orc <- .nlp_oracles(prob, use_hessian = FALSE)
  hs <- orc$hessianstructure()
  expect_equal(length(hs$rows), 0L)
  expect_error(orc$hessian(c(1, 1), duals = numeric(0), obj_factor = 1), "use_hessian")
})

# ===================================================================
# End-to-end nlp solve -- pending the Uno backend (Step 6)
# ===================================================================

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_analytic_polytope_center
test_that("nlp=TRUE solve via Uno", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## A DNLP whose data the diff_engine + Uno handle end-to-end:
  ## min sum(exp(x))  s.t.  sum(x) == 1.5  =>  x_i = 0.5, f* = 3 exp(0.5).
  x <- Variable(3)
  value(x) <- c(0, 0.5, 1)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
  val <- psolve(prob, solver = "UNO", nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 3 * exp(0.5), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(0.5, 3), tolerance = 1e-3)
})
