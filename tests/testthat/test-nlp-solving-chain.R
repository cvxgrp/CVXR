# Tests for the NLP solving chain (CVXPY 1.9.0, Wave-1 Step 6).
#
# nlp_solving_chain.R (.build_nlp_chain / .set_nlp_initial_point / solve_nlp),
# defines.R (SOLVER_MAP_NLP / NLP_SOLVER_VARIANTS / .nlp_solver_available),
# uno_nlpif.R (UNO(NLPsolver)), and the psolve(nlp = TRUE) dispatch.
#
# CVXPY SOURCE: reductions/solvers/nlp_solving_chain.py,
#   reductions/solvers/defines.py, reductions/solvers/nlp_solvers/uno_nlpif.py,
#   problems/problem.py:1109-1115

library(CVXR)

# ===================================================================
# Solver registry + variants (no solver invocation needed)
# ===================================================================

## @cvxpy defines.py
test_that("SOLVER_MAP_NLP holds the four NLP solvers in preference order", {
  expect_equal(NLP_SOLVER_PREFERENCE, c("IPOPT", "KNITRO", "UNO", "COPT"))
  ## IPOPT and UNO are available when their optional R packages are installed.
  expect_equal(.nlp_solver_available("IPOPT"),
               requireNamespace("ipopt", quietly = TRUE))
  expect_equal(.nlp_solver_available("UNO"),
               requireNamespace("Uno", quietly = TRUE))
  expect_equal(.nlp_solver_available("KNITRO"), FALSE)
  expect_equal(.nlp_solver_available("COPT"), FALSE)
})

## @cvxpy defines.py
test_that("NLP_SOLVER_VARIANTS maps variant names to base + kwargs", {
  expect_equal(NLP_SOLVER_VARIANTS$uno_sqp$base, "UNO")
  expect_equal(NLP_SOLVER_VARIANTS$uno_sqp$kwargs$preset, "filtersqp")
  expect_equal(NLP_SOLVER_VARIANTS$uno_ipm$kwargs$preset, "ipopt")
  expect_equal(NLP_SOLVER_VARIANTS$uno_ipm$kwargs$linear_solver, "MUMPS")
  expect_equal(NLP_SOLVER_VARIANTS$knitro_sqp$base, "KNITRO")
})

# ===================================================================
# psolve(nlp = TRUE) dispatch + error paths
# ===================================================================

## @cvxpy problem.py:1114-1115 (DNLPError)
test_that("nlp = TRUE on a non-DNLP problem errors", {
  ## A difference of piecewise-linear convex atoms is neither DCP nor smooth,
  ## so it is not DNLP (unlike a difference of smooth atoms such as
  ## sum_squares(x) - sum_squares(y), which is smooth and therefore DNLP).
  x <- Variable(2); y <- Variable(2)
  prob <- Problem(Minimize(sum_entries(abs(x)) - sum_entries(abs(y))))
  expect_false(is_dnlp(prob))
  expect_error(psolve(prob, nlp = TRUE), "not DNLP")
})

## @cvxpy NONE
test_that("an explicit unavailable NLP solver errors", {
  skip_if_not_installed("sparsediff")
  x <- Variable(2)
  value(x) <- c(0, 0)
  prob <- Problem(Minimize(sum_squares(x - 1)), list(sum_entries(x) <= 10))
  expect_error(psolve(prob, solver = "KNITRO", nlp = TRUE), "not available")
})

# ===================================================================
# End-to-end DNLP solves via Uno
# ===================================================================

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_qcp
test_that("nlp = TRUE solves a convex QP-shaped DNLP (auto-select NLP solver)", {
  skip_if_not_installed("sparsediff")
  skip_if_not(.nlp_solver_available("IPOPT") || .nlp_solver_available("UNO"),
              "no optional NLP solver installed")
  ## min ||x - (1,2)||^2  s.t.  sum(x) <= 10.  Unconstrained optimum (1,2) is
  ## feasible, so x* = (1,2), f* = 0.  DCP => DNLP.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(sum_entries(x) <= 10))
  val <- psolve(prob, nlp = TRUE, print_level = 0L)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-4)
})

## @cvxpy reductions/solvers/defines.py reductions/solvers/nlp_solving_chain.py
test_that("nlp = TRUE accepts explicit IPOPT", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("ipopt")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(sum_entries(x) <= 10))
  expect_equal(psolve(prob, solver = "IPOPT", nlp = TRUE, print_level = 0L),
               0, tolerance = 1e-5)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-4)
})

## @cvxpy reductions/solvers/defines.py reductions/solvers/nlp_solving_chain.py
test_that("nlp = TRUE accepts explicit UNO and a variant name", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(sum_entries(x) <= 10))
  expect_equal(psolve(prob, solver = "UNO", nlp = TRUE), 0, tolerance = 1e-5)

  ## A variant name (uno_sqp) resolves to the UNO base with its preset kwargs.
  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_squares(x2 - c(1, 2))),
                   list(sum_entries(x2) <= 10))
  expect_equal(psolve(prob2, solver = "uno_sqp", nlp = TRUE), 0, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_hyperbolic.py::TestHyperbolic::test_sinh
test_that("uno_ipm variant solves an indefinite-Hessian DNLP", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## min sum(sinh(logistic(x * 2)))  s.t.  x >= 0.1, sum(x) == 10.  The filtersqp
  ## preset's QP subproblem solver (HiGHS in this build) cannot factor the
  ## indefinite Lagrangian Hessian here ("Algorithmic error"), so this exercises
  ## the "uno_ipm" variant (interior point + MUMPS) -- the same engine CVXR's UNO
  ## now defaults to -- which converges to the optimum (x_i = 1).
  n <- 10
  x <- Variable(n)
  prob <- Problem(Minimize(sum_entries(sinh(logistic(x * 2)))),
                  list(x >= 0.1, sum_entries(x) == 10))
  val <- psolve(prob, solver = "uno_ipm", nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 41.349266, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_least_squares
test_that("nlp = TRUE solves a problem with variable bounds (sign attribute)", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## min sum(exp(x))  s.t.  sum(x) == 1.5,  x >= 0 (nonneg).  By symmetry the
  ## optimum is x_i = 0.5, f* = 3 * exp(0.5).
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(exp(x))), list(sum_entries(x) == 1.5))
  val <- psolve(prob, nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 3 * exp(0.5), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(0.5, 3), tolerance = 1e-3)
})
