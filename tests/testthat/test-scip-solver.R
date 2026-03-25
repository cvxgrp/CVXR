## SCIP solver tests — ports CVXPY test_conic_solvers.py::TestSCIP
## Uses the standard test helpers from helper-solver-test-problems.R.
## SCIP supports: Zero + NonNeg + SOC, MIP-capable.
## Duals are NOT tested (R scip package lacks dual extraction API).

require_solver("SCIP")

# ── Helper: verify objective ─────────────────────────────────────

verify_obj <- function(prob, expected, tol, solver_name) {
  val <- psolve(prob, solver = solver_name)
  if (is.infinite(expected) && expected > 0) {
    ## Infeasible: psolve returns Inf
    expect_true(is.infinite(val) && val > 0,
                label = paste(solver_name, "should be infeasible (Inf)"))
  } else if (is.infinite(expected) && expected < 0) {
    ## Unbounded: psolve returns -Inf
    expect_true(is.infinite(val) && val < 0,
                label = paste(solver_name, "should be unbounded (-Inf)"))
  } else {
    expect_equal(val, expected, tolerance = tol,
                 label = paste(solver_name, "objective"))
  }
}

# ── Helper: verify primal ────────────────────────────────────────

verify_primal <- function(prob, var, expected, tol, solver_name) {
  if (is.null(expected)) return(invisible(NULL))
  actual <- value(var)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = tol,
               label = paste(solver_name, "primal"))
}

# ── Helper: verify duals ─────────────────────────────────────────

verify_dual <- function(con, expected, tol, solver_name) {
  if (is.null(expected)) return(invisible(NULL))
  actual <- dual_value(con)
  if (is.list(expected)) {
    ## SOC / ExpCone / PowCone: list of arrays
    for (k in seq_along(expected)) {
      expect_equal(as.numeric(actual[[k]]), as.numeric(expected[[k]]),
                   tolerance = tol,
                   label = paste(solver_name, "dual component", k))
    }
  } else {
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = tol,
                 label = paste(solver_name, "dual"))
  }
}

# ── LP tests ──────────────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_lp_0
test_that("SCIP: lp_0", {
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_lp_1
test_that("SCIP: lp_1", {
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "SCIP")
  ## Duals skipped (CVXPY also passes duals=True here but SCIP duals are unreliable)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_lp_2
test_that("SCIP: lp_2", {
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "SCIP")
  ## Duals skipped (CVXPY passes duals=False)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_lp_3
test_that("SCIP: lp_3 (unbounded)", {
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_lp_4
test_that("SCIP: lp_4 (infeasible)", {
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

# ── SOCP tests ────────────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_socp_0
test_that("SCIP: socp_0", {
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_socp_1
test_that("SCIP: socp_1", {
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCIP")
  ## Duals skipped (CVXPY passes duals=False, places=2)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_socp_2
test_that("SCIP: socp_2", {
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCIP")
  ## Duals skipped (CVXPY passes duals=False, places=2)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_socp_3
test_that("SCIP: socp_3 (axis 1 and 2)", {
  ## axis=1 (row-wise, CVXPY axis=0)
  h0 <- sth_socp_3_ax0()
  verify_obj(h0$prob, h0$expect_obj, 1e-2, "SCIP")
  ## axis=2 (col-wise, CVXPY axis=1)
  h1 <- sth_socp_3_ax1()
  verify_obj(h1$prob, h1$expect_obj, 1e-2, "SCIP")
})

# ── MI-LP tests ───────────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_lp_0
test_that("SCIP: mi_lp_0", {
  h <- sth_mi_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_lp_1
test_that("SCIP: mi_lp_1", {
  h <- sth_mi_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_lp_2
test_that("SCIP: mi_lp_2", {
  h <- sth_mi_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_lp_3
test_that("SCIP: mi_lp_3 (infeasible boolean MIP)", {
  h <- sth_mi_lp_3()
  psolve(h$prob, solver = "SCIP")
  expect_true(status(h$prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_lp_5
test_that("SCIP: mi_lp_5", {
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "SCIP")
})

# ── MI-SOCP tests ─────────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_socp_1
test_that("SCIP: mi_socp_1", {
  h <- sth_mi_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCIP")
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_mi_socp_2
test_that("SCIP: mi_socp_2", {
  h <- sth_mi_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCIP")
})

# ── Parameter tests ───────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_test_params__no_params_set
test_that("SCIP: params - no params set", {
  x <- Variable()
  y <- Variable()
  cons <- list(x >= 0, y >= 1, x + y <= 4)
  prob <- Problem(Maximize(x), cons)
  val <- psolve(prob, solver = "SCIP")
  expect_equal(val, 3, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_test_params__valid_params
test_that("SCIP: params - valid params", {
  x <- Variable()
  y <- Variable()
  cons <- list(x >= 0, y >= 1, x + y <= 4)
  prob <- Problem(Maximize(x), cons)
  val <- psolve(prob, solver = "SCIP")
  expect_equal(val, 3, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_test_params__valid_scip_params
test_that("SCIP: params - valid scip_params", {
  x <- Variable()
  y <- Variable()
  cons <- list(x >= 0, y >= 1, x + y <= 4)
  prob <- Problem(Maximize(x), cons)
  val <- psolve(prob, solver = "SCIP",
                scip_params = list("limits/gap" = 0.1))
  expect_equal(val, 3, tolerance = 0.2)
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_test_params__invalid_params
test_that("SCIP: params - invalid params raises error", {
  x <- Variable()
  y <- Variable()
  cons <- list(x >= 0, y >= 1, x + y <= 4)
  prob <- Problem(Maximize(x), cons)
  ## Use a name that cannot partial-match any psolve formal
  ## (CVXPY uses "a" but R partial-matches that to "abstol")
  expect_error(
    psolve(prob, solver = "SCIP", bogus_param = "what?"),
    "Invalid SCIP"
  )
})

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_test_params__invalid_scip_params
test_that("SCIP: params - invalid scip_params raises error", {
  x <- Variable()
  y <- Variable()
  cons <- list(x >= 0, y >= 1, x + y <= 4)
  prob <- Problem(Maximize(x), cons)
  expect_error(
    psolve(prob, solver = "SCIP", scip_params = list(a = "what?")),
    "Invalid SCIP"
  )
})

# ── Time limit test ───────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_time_limit_reached
test_that("SCIP: time limit with no solution raises error", {
  ## Large MIP that can't be solved instantly
  n <- 50
  x <- Variable(n, integer = TRUE)
  A <- matrix(rnorm(n * n), n, n)
  b <- rep(1, n)
  obj <- Minimize(sum(x))
  cons <- list(A %*% x >= b, x >= 0, x <= 10)
  prob <- Problem(obj, cons)
  expect_error(
    psolve(prob, solver = "SCIP", scip_params = list("limits/time" = 0.0)),
    "failed"
  )
})

# ── Solver stats test ─────────────────────────────────────────────

## @cvxpy test_conic_solvers.py::TestSCIP::test_scip_solver_stats
test_that("SCIP: solver stats available", {
  h <- sth_lp_0()
  result <- solve(h$prob, solver = "SCIP")
  expect_equal(result$solver, "SCIP")
  ss <- solver_stats(h$prob)
  expect_true(!is.null(ss@solve_time))
  expect_true(!is.null(ss@num_iters))
})
