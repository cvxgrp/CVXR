## Tests for standard solver parameter translation
## Tests that solver_default_param() and .apply_std_params() work correctly,
## and that psolve() properly forwards standard params to solvers.

context("Standard Solver Parameters")

# ── Unit tests for .apply_std_params() ─────────────────────────────

## @cvxpy NONE
test_that(".apply_std_params translates standard names to Clarabel native names", {
  opts <- CVXR:::.apply_std_params("CLARABEL", list(),
                                   feastol = 1e-3, reltol = 1e-4,
                                   abstol = 1e-5, num_iter = 100L)
  expect_equal(opts$tol_feas, 1e-3)
  expect_equal(opts$tol_gap_rel, 1e-4)
  expect_equal(opts$tol_gap_abs, 1e-5)
  expect_equal(opts$max_iter, 100L)
})

## @cvxpy NONE
test_that(".apply_std_params translates standard names to OSQP native names", {
  opts <- CVXR:::.apply_std_params("OSQP", list(),
                                   feastol = 1e-3, reltol = 1e-4,
                                   abstol = 1e-5, num_iter = 500L)
  expect_equal(opts$eps_prim_inf, 1e-3)
  expect_equal(opts$eps_rel, 1e-4)
  expect_equal(opts$eps_abs, 1e-5)
  expect_equal(opts$max_iter, 500L)
})

## @cvxpy NONE
test_that(".apply_std_params translates standard names to SCS native names", {
  opts <- CVXR:::.apply_std_params("SCS", list(),
                                   feastol = 1e-3, reltol = 1e-4,
                                   abstol = 1e-5, num_iter = 1000L)
  ## SCS has no feastol mapping — should be absent

  expect_null(opts$feastol)
  expect_equal(opts$eps_rel, 1e-4)
  expect_equal(opts$eps_abs, 1e-5)
  expect_equal(opts$max_iters, 1000L)
})

## @cvxpy NONE
test_that(".apply_std_params translates standard names to ECOS native names", {
  opts <- CVXR:::.apply_std_params("ECOS", list(),
                                   feastol = 1e-6, reltol = 1e-7,
                                   abstol = 1e-7, num_iter = 200L)
  expect_equal(opts$FEASTOL, 1e-6)
  expect_equal(opts$RELTOL, 1e-7)
  expect_equal(opts$ABSTOL, 1e-7)
  expect_equal(opts$MAXIT, 200L)
})

## @cvxpy NONE
test_that(".apply_std_params does not override solver-native names in opts", {
  ## User passes both feastol=1e-3 (standard) and tol_feas=1e-6 (native)
  ## Native should win
  opts <- CVXR:::.apply_std_params("CLARABEL",
                                   list(tol_feas = 1e-6),
                                   feastol = 1e-3, reltol = NULL,
                                   abstol = NULL, num_iter = NULL)
  expect_equal(opts$tol_feas, 1e-6)
})

## @cvxpy NONE
test_that(".apply_std_params leaves opts unchanged for NULL standard params", {
  opts <- CVXR:::.apply_std_params("CLARABEL", list(extra = 42),
                                   feastol = NULL, reltol = NULL,
                                   abstol = NULL, num_iter = NULL)
  expect_equal(opts, list(extra = 42))
})

## @cvxpy NONE
test_that(".apply_std_params returns opts unchanged for unknown solver", {
  opts <- CVXR:::.apply_std_params("UNKNOWN_SOLVER", list(foo = 1),
                                   feastol = 1e-3, reltol = 1e-4,
                                   abstol = 1e-5, num_iter = 100L)
  expect_equal(opts, list(foo = 1))
})

## @cvxpy NONE
test_that(".apply_std_params returns opts unchanged for GLPK (no mappings)", {
  opts <- CVXR:::.apply_std_params("GLPK", list(),
                                   feastol = 1e-3, reltol = 1e-4,
                                   abstol = 1e-5, num_iter = 100L)
  expect_equal(opts, list())
})

# ── Integration tests: psolve() with standard params ──────────────

## @cvxpy NONE
test_that("psolve with feastol works for Clarabel", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  result <- psolve(prob, solver = "CLARABEL", feastol = 1e-3)
  expect_equal(result, 2.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("psolve with num_iter works for Clarabel", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  result <- psolve(prob, solver = "CLARABEL", num_iter = 500)
  expect_equal(result, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("psolve with reltol and abstol works for SCS", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  result <- psolve(prob, solver = "SCS", reltol = 1e-6, abstol = 1e-6)
  expect_equal(result, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("psolve with num_iter works for OSQP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  result <- psolve(prob, solver = "OSQP", num_iter = 20000)
  expect_equal(result, 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("psolve with num_iter works for ECOS", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  result <- psolve(prob, solver = "ECOS", num_iter = 200)
  expect_equal(result, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("solver-native params in ... override standard params", {
  ## Pass both feastol (standard) and tol_feas (native); native should win
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  ## This should use tol_feas=1e-8 (native), not feastol=1e-1 (standard)
  result <- psolve(prob, solver = "CLARABEL", feastol = 1e-1, tol_feas = 1e-8)
  expect_equal(result, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("solver_default_param() is exported and has expected structure", {
  expect_true(is.list(solver_default_param()))
  expect_true("CLARABEL" %in% names(solver_default_param()))
  expect_true("OSQP" %in% names(solver_default_param()))
  expect_true("SCS" %in% names(solver_default_param()))
  expect_true("ECOS" %in% names(solver_default_param()))

  ## Check Clarabel structure
  cl <- solver_default_param()$CLARABEL
  expect_true(all(c("reltol", "abstol", "feastol", "num_iter") %in% names(cl)))
  expect_equal(cl$feastol$name, "tol_feas")
  expect_equal(cl$feastol$value, 1e-8)
})
