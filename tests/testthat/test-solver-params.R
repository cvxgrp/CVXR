## Tests for standard solver parameter translation
## Tests that solver_default_param() and .build_solver_params() work
## correctly, and that psolve() properly forwards standard params to
## solvers. (These unit tests originally targeted .apply_std_params(),
## a duplicate of the translation rule with no production callers; it
## was deleted and the tests re-pointed at the live path.)

context("Standard Solver Parameters")

# ── Unit tests for .build_solver_params() ──────────────────────────

## @cvxpy NONE
test_that("standard names translate to Clarabel native names", {
  opts <- CVXR:::.build_solver_params("CLARABEL",
                                      solver_opts(feastol = 1e-3, reltol = 1e-4,
                                                  abstol = 1e-5, num_iter = 100L))
  expect_equal(opts$tol_feas, 1e-3)
  expect_equal(opts$tol_gap_rel, 1e-4)
  expect_equal(opts$tol_gap_abs, 1e-5)
  expect_equal(opts$max_iter, 100L)
})

## @cvxpy NONE
test_that("standard names translate to OSQP native names", {
  opts <- CVXR:::.build_solver_params("OSQP",
                                      solver_opts(feastol = 1e-3, reltol = 1e-4,
                                                  abstol = 1e-5, num_iter = 500L))
  expect_equal(opts$eps_prim_inf, 1e-3)
  expect_equal(opts$eps_rel, 1e-4)
  expect_equal(opts$eps_abs, 1e-5)
  expect_equal(opts$max_iter, 500L)
})

## @cvxpy NONE
test_that("standard names translate to SCS native names", {
  ## SCS has no feastol mapping — dropped, and since the unmapped-param
  ## warning landed, dropped LOUDLY.
  expect_warning(
    opts <- CVXR:::.build_solver_params("SCS",
                                        solver_opts(feastol = 1e-3, reltol = 1e-4,
                                                    abstol = 1e-5, num_iter = 1000L)),
    class = "cvxr_unmapped_solver_param"
  )
  expect_null(opts$feastol)
  expect_equal(opts$eps_rel, 1e-4)
  expect_equal(opts$eps_abs, 1e-5)
  expect_equal(opts$max_iters, 1000L)
})

## @cvxpy NONE
test_that("standard names translate to ECOS native names", {
  opts <- CVXR:::.build_solver_params("ECOS",
                                      solver_opts(feastol = 1e-6, reltol = 1e-7,
                                                  abstol = 1e-7, num_iter = 200L))
  expect_equal(opts$FEASTOL, 1e-6)
  expect_equal(opts$RELTOL, 1e-7)
  expect_equal(opts$ABSTOL, 1e-7)
  expect_equal(opts$MAXIT, 200L)
})

## @cvxpy NONE
test_that("translation does not override solver-native names in opts", {
  ## User passes both feastol=1e-3 (standard) and tol_feas=1e-6 (native)
  ## Native should win
  opts <- CVXR:::.build_solver_params("CLARABEL",
                                      solver_opts(feastol = 1e-3, tol_feas = 1e-6))
  expect_equal(opts$tol_feas, 1e-6)
})

## @cvxpy NONE
test_that("opts pass through unchanged when no standard params are set", {
  opts <- CVXR:::.build_solver_params("CLARABEL", solver_opts(extra = 42))
  expect_equal(opts, list(extra = 42))
})

## @cvxpy NONE
test_that("an unknown solver drops standard params, with a warning", {
  expect_warning(
    opts <- CVXR:::.build_solver_params("UNKNOWN_SOLVER",
                                        solver_opts(feastol = 1e-3, reltol = 1e-4,
                                                    abstol = 1e-5, num_iter = 100L,
                                                    foo = 1)),
    class = "cvxr_unmapped_solver_param"
  )
  expect_equal(opts, list(foo = 1))
})

## @cvxpy NONE
test_that("GLPK (no mappings) drops standard params, with a warning", {
  expect_warning(
    opts <- CVXR:::.build_solver_params("GLPK",
                                        solver_opts(feastol = 1e-3, reltol = 1e-4,
                                                    abstol = 1e-5, num_iter = 100L)),
    class = "cvxr_unmapped_solver_param"
  )
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

## @cvxpy NONE
test_that(".apply_std_params is gone -- one translation path only", {
  expect_false(exists(".apply_std_params", envir = asNamespace("CVXR")))
})
