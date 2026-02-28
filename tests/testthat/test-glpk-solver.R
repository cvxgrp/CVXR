## Tests for GLPK and GLPK_MI solver integration
## All tests are wrapped in skip_if_not_installed("Rglpk") so they
## pass cleanly when Rglpk is not available.

# ══════════════════════════════════════════════════════════════════
# Unit tests — constants and availability
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK_SOLVER constant is exported", {
  expect_equal(CVXR::GLPK_SOLVER, "GLPK")
})

## @cvxpy NONE
test_that("GLPK_MI_SOLVER constant is exported", {
  expect_equal(CVXR::GLPK_MI_SOLVER, "GLPK_MI")
})

## @cvxpy NONE
test_that("GLPK appears in installed_solvers when available", {
  skip_if_not_installed("Rglpk")
  solvers <- installed_solvers()
  expect_true("GLPK" %in% solvers)
  expect_true("GLPK_MI" %in% solvers)
})

# ══════════════════════════════════════════════════════════════════
# GLPK LP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK basic LP", {
  skip_if_not_installed("Rglpk")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GLPK")
  expect_equal(val, 3.0, tolerance = 1e-6)
  expect_equal(value(x), matrix(c(3.0, 0.0), ncol = 1), tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK LP with equality constraint", {
  skip_if_not_installed("Rglpk")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] == 5, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "GLPK")
  expect_equal(val, 5.0, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK Maximize", {
  skip_if_not_installed("Rglpk")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GLPK")
  expect_equal(val, 10.0, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK infeasible LP", {
  skip_if_not_installed("Rglpk")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x >= 5, x <= 1))
  psolve(prob, solver = "GLPK")
  expect_equal(status(prob), "infeasible")
})

## @cvxpy NONE
test_that("GLPK unbounded LP", {
  skip_if_not_installed("Rglpk")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "GLPK")
  expect_equal(status(prob), "unbounded")
})

# ══════════════════════════════════════════════════════════════════
# GLPK dual values
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK LP duals match Clarabel", {
  skip_if_not_installed("Rglpk")
  skip_if_not_installed("clarabel")

  ## Solve with GLPK (fresh objects — no cache leakage)
  x1 <- Variable(name = "x1")
  y1 <- Variable(name = "y1")
  eq1 <- (x1 + y1 == 3)
  ineq1a <- (x1 >= 1)
  ineq1b <- (y1 >= 0)
  prob1 <- Problem(Minimize(x1 + 2 * y1), list(eq1, ineq1a, ineq1b))
  psolve(prob1, solver = "GLPK")

  ## Solve with Clarabel (fresh objects)
  x2 <- Variable(name = "x2")
  y2 <- Variable(name = "y2")
  eq2 <- (x2 + y2 == 3)
  ineq2a <- (x2 >= 1)
  ineq2b <- (y2 >= 0)
  prob2 <- Problem(Minimize(x2 + 2 * y2), list(eq2, ineq2a, ineq2b))
  psolve(prob2, solver = "CLARABEL")

  ## Compare duals
  expect_equal(dual_value(eq1), dual_value(eq2), tolerance = 1e-4)
  expect_equal(dual_value(ineq1a), dual_value(ineq2a), tolerance = 1e-4)
  expect_equal(dual_value(ineq1b), dual_value(ineq2b), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("GLPK LP duals verified against CVXPY reference", {
  skip_if_not_installed("Rglpk")
  ## CVXPY GLPK reference (via uv run python):
  ##   eq dual = -1.0, ineq1 dual ≈ 0, ineq2 dual = 1.0
  ## Our conic path duals should match Clarabel/SCS convention:
  ##   eq dual = -1, ineq (x>=1, inactive) ≈ 0, ineq (y>=0, active) = 1
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  eq <- (x + y == 3)
  ineq1 <- (x >= 1)
  ineq2 <- (y >= 0)
  prob <- Problem(Minimize(x + 2 * y), list(eq, ineq1, ineq2))
  psolve(prob, solver = "GLPK")

  expect_equal(value(prob), 3.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(y)), 0.0, tolerance = 1e-6)
  expect_equal(as.numeric(dual_value(eq)), -1.0, tolerance = 1e-6)
  expect_equal(as.numeric(dual_value(ineq1)), 0.0, tolerance = 1e-6)
  expect_equal(as.numeric(dual_value(ineq2)), 1.0, tolerance = 1e-6)
})

# ══════════════════════════════════════════════════════════════════
# GLPK_MI MILP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK_MI basic integer LP", {
  skip_if_not_installed("Rglpk")
  x <- Variable(integer = TRUE, name = "x_int")
  y <- Variable(name = "y_cont")
  prob <- Problem(Minimize(x + 2 * y),
                  list(x + y == 3, x >= 1, y >= 0))
  val <- psolve(prob, solver = "GLPK_MI")
  expect_equal(val, 3.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(y)), 0.0, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK_MI boolean variable", {
  skip_if_not_installed("Rglpk")
  b <- Variable(boolean = TRUE, name = "b")
  x <- Variable(name = "x_cont")
  prob <- Problem(Minimize(x), list(x >= b, x >= 0.5))
  val <- psolve(prob, solver = "GLPK_MI")
  expect_equal(val, 0.5, tolerance = 1e-6)
  expect_equal(as.numeric(value(b)), 0.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK_MI integer LP where integrality matters", {
  skip_if_not_installed("Rglpk")
  ## min -x s.t. x <= 3.5, x integer, x >= 0
  ## Continuous optimum: x = 3.5, integer optimum: x = 3
  x <- Variable(integer = TRUE, name = "x_int")
  prob <- Problem(Maximize(x), list(x <= 3.5, x >= 0))
  val <- psolve(prob, solver = "GLPK_MI")
  expect_equal(val, 3.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-6)
})

## @cvxpy NONE
test_that("GLPK_MI does not return duals", {
  skip_if_not_installed("Rglpk")
  x <- Variable(integer = TRUE, name = "x_int")
  constr <- (x >= 1)
  prob <- Problem(Minimize(x), list(constr))
  psolve(prob, solver = "GLPK_MI")
  ## MIP problems should not have dual values
  expect_null(dual_value(constr))
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK matches SCS on LP", {
  skip_if_not_installed("Rglpk")
  skip_if_not_installed("scs")

  ## GLPK
  x1 <- Variable(3)
  prob1 <- Problem(Minimize(x1[1] + x1[2] + x1[3]),
                   list(x1 >= 1, x1[1] + x1[2] <= 5))
  val1 <- psolve(prob1, solver = "GLPK")

  ## SCS
  x2 <- Variable(3)
  prob2 <- Problem(Minimize(x2[1] + x2[2] + x2[3]),
                   list(x2 >= 1, x2[1] + x2[2] <= 5))
  val2 <- psolve(prob2, solver = "SCS")

  expect_equal(val1, val2, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("GLPK and GLPK_MI give same result on pure LP", {
  skip_if_not_installed("Rglpk")

  ## GLPK (LP solver)
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(x1[1] + 3 * x1[2]),
                   list(x1[1] + x1[2] >= 4, x1 >= 0))
  val1 <- psolve(prob1, solver = "GLPK")

  ## GLPK_MI on same continuous LP
  x2 <- Variable(2)
  prob2 <- Problem(Minimize(x2[1] + 3 * x2[2]),
                   list(x2[1] + x2[2] >= 4, x2 >= 0))
  val2 <- psolve(prob2, solver = "GLPK_MI")

  expect_equal(val1, val2, tolerance = 1e-6)
})

# ══════════════════════════════════════════════════════════════════
# Larger LP problem
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GLPK larger LP problem", {
  skip_if_not_installed("Rglpk")
  n <- 20
  set.seed(42)
  c_vec <- rnorm(n)
  A_mat <- matrix(rnorm(n * n), nrow = n)

  x <- Variable(n)
  prob <- Problem(Minimize(t(c_vec) %*% x),
                  list(A_mat %*% x <= rep(10, n), x >= -5, x <= 5))
  val <- psolve(prob, solver = "GLPK")
  expect_equal(status(prob), "optimal")

  ## Cross-check with Clarabel
  skip_if_not_installed("clarabel")
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(t(c_vec) %*% x2),
                   list(A_mat %*% x2 <= rep(10, n), x2 >= -5, x2 <= 5))
  val2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val, val2, tolerance = 1e-4)
})
