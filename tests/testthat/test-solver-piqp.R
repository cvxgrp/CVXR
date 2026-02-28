## Tests for PIQP solver integration (backed by R piqp package)
## All tests are wrapped in skip_if_not_installed("piqp") so they
## pass cleanly when piqp is not available.

# ══════════════════════════════════════════════════════════════════
# Unit tests — constants and availability
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP_SOLVER constant is exported", {
  expect_equal(CVXR::PIQP_SOLVER, "PIQP")
})

## @cvxpy NONE
test_that("PIQP appears in installed_solvers when available", {
  skip_if_not_installed("piqp")
  solvers <- installed_solvers()
  expect_true("PIQP" %in% solvers)
})

# ══════════════════════════════════════════════════════════════════
# PIQP LP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP basic LP", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 1, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "PIQP")
  ## CVXPY ref: x = [1, 0], obj = 1
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1.0, 0.0), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("PIQP LP with equality constraint", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] == 5, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "PIQP")
  expect_equal(val, 5.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("PIQP Maximize LP", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "PIQP")
  expect_equal(val, 10.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("PIQP lp_1 (classic example with duals)", {
  ## min -4x₀ - 5x₁, s.t. 2x₀ + x₁ ≤ 3, x₀ + 2x₁ ≤ 3, x ≥ 0
  ## Reference: x = [1, 1], obj = -9, duals = [1, 2, 0, 0]
  skip_if_not_installed("piqp")
  x <- Variable(2, name = "x")
  c1 <- (2 * x[1] + x[2] <= 3)
  c2 <- (x[1] + 2 * x[2] <= 3)
  c3 <- (x[1] >= 0)
  c4 <- (x[2] >= 0)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]), list(c1, c2, c3, c4))
  psolve(prob, solver = "PIQP")
  expect_equal(value(prob), -9.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), 0.0, tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(c4)), 0.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("PIQP lp_2 (bounds LP with duals)", {
  ## min x₀ + 0.5x₁, s.t. -100 ≤ x₀ ≤ -10, x₁ == 1
  ## Reference: x = [-100, 1], obj = -99.5, duals = [1, 0, -0.5]
  skip_if_not_installed("piqp")
  x <- Variable(2, name = "x")
  c1 <- (x[1] >= -100)
  c2 <- (x[1] <= -10)
  c3 <- (x[2] == 1)
  prob <- Problem(Minimize(x[1] + 0.5 * x[2]), list(c1, c2, c3))
  psolve(prob, solver = "PIQP")
  expect_equal(value(prob), -99.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(-100, 1), tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), -0.5, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# PIQP QP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP simple QP", {
  skip_if_not_installed("piqp")
  ## CVXPY ref: min quad_form(x, P) + x₁ + x₂, s.t. x₁ + x₂ == 1, x ≥ 0
  ## P = [[2, 0.5], [0.5, 1]], x = [0.25, 0.75], obj = 1.875
  x <- Variable(2)
  P_data <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P_data) + x[1] + x[2]),
                  list(x[1] + x[2] == 1, x >= 0))
  val <- psolve(prob, solver = "PIQP")
  expect_equal(val, 1.875, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0.25, 0.75), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("PIQP QP with duals (sum_squares)", {
  ## min sum_squares(x) s.t. sum(x) == 1, x >= 0
  ## Reference: x = [1/3, 1/3, 1/3], obj = 1/3, eq_dual = -2/3
  skip_if_not_installed("piqp")
  x <- Variable(3)
  eq <- (sum_entries(x) == 1)
  ineq <- (x >= 0)
  prob <- Problem(Minimize(sum_squares(x)), list(eq, ineq))
  psolve(prob, solver = "PIQP")
  expect_equal(value(prob), 1/3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(eq)), -2/3, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ineq)), rep(0, 3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("PIQP QP (quadratic objective, no constraints)", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(3, 4))))
  val <- psolve(prob, solver = "PIQP")
  expect_equal(val, 0.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(3, 4), tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# PIQP infeasible / unbounded
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP infeasible problem", {
  skip_if_not_installed("piqp")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x >= 5, x <= 1))
  psolve(prob, solver = "PIQP")
  ## PIQP may return infeasible or user_limit for infeasible problems
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate",
                                    "user_limit", "solver_error"))
})

## @cvxpy NONE
test_that("PIQP unbounded problem", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 1))
  psolve(prob, solver = "PIQP")
  ## PIQP may return unbounded or user_limit for unbounded problems
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate",
                                    "user_limit", "solver_error"))
})

# ══════════════════════════════════════════════════════════════════
# PIQP dual values (LP)
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP dual values (LP)", {
  skip_if_not_installed("piqp")
  ## min x + 2y  s.t.  x + y == 3, x >= 0, y >= 0
  ## Reference: value = 3, x = [3, 0]
  ## eq dual = -1, ineq duals = [~0, 1]
  x <- Variable(2)
  eq <- (x[1] + x[2] == 3)
  ineq <- (x >= 0)
  prob <- Problem(Minimize(x[1] + 2 * x[2]), list(eq, ineq))
  psolve(prob, solver = "PIQP")

  expect_equal(value(prob), 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(eq)), -1.0, tolerance = 1e-3)
  d <- as.numeric(dual_value(ineq))
  expect_equal(d[1], 0.0, tolerance = 1e-3)
  expect_equal(d[2], 1.0, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# PIQP solver options + standard params
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP solver_opts (tolerances)", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "PIQP",
                abstol = 1e-7, reltol = 1e-7, num_iter = 500L)
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("PIQP native solver_opts", {
  skip_if_not_installed("piqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "PIQP",
                max_iter = 500L, eps_abs = 1e-7)
  expect_equal(val, 1.0, tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# PIQP warm-start
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP warm-start", {
  skip_if_not_installed("piqp")
  ## Solve a parametric QP twice, second solve should warm-start
  a <- Parameter()
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + a * x[1]),
                  list(x >= 0, sum_entries(x) == 1))

  ## First solve
  value(a) <- 1
  val1 <- psolve(prob, solver = "PIQP", warm_start = TRUE)
  x1 <- as.numeric(value(x))

  ## Second solve (warm-started)
  value(a) <- 2
  val2 <- psolve(prob, solver = "PIQP", warm_start = TRUE)
  x2 <- as.numeric(value(x))

  expect_true(val2 > val1)
  expect_true(x2[2] > x1[2])
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity: PIQP vs OSQP
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PIQP vs OSQP: LP value + primals + duals", {
  skip_if_not_installed("piqp")
  skip_if_not_installed("osqp")

  ## PIQP (fresh objects)
  x1 <- Variable(2)
  eq1 <- (x1[1] + x1[2] == 3)
  ge1a <- (x1[1] >= 1)
  ge1b <- (x1[2] >= 0)
  prob1 <- Problem(Minimize(x1[1] + 2 * x1[2]), list(eq1, ge1a, ge1b))
  val1 <- psolve(prob1, solver = "PIQP")

  ## OSQP (fresh objects)
  x2 <- Variable(2)
  eq2 <- (x2[1] + x2[2] == 3)
  ge2a <- (x2[1] >= 1)
  ge2b <- (x2[2] >= 0)
  prob2 <- Problem(Minimize(x2[1] + 2 * x2[2]), list(eq2, ge2a, ge2b))
  val2 <- psolve(prob2, solver = "OSQP")

  expect_equal(val1, val2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(eq1)), as.numeric(dual_value(eq2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1a)), as.numeric(dual_value(ge2a)), tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(ge1b)), as.numeric(dual_value(ge2b)), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("PIQP vs OSQP: QP value + primals", {
  skip_if_not_installed("piqp")
  skip_if_not_installed("osqp")

  ## PIQP
  x1 <- Variable(3)
  prob1 <- Problem(Minimize(sum_squares(x1)), list(sum_entries(x1) == 1, x1 >= 0))
  val1 <- psolve(prob1, solver = "PIQP")

  ## OSQP
  x2 <- Variable(3)
  prob2 <- Problem(Minimize(sum_squares(x2)), list(sum_entries(x2) == 1, x2 >= 0))
  val2 <- psolve(prob2, solver = "OSQP")

  expect_equal(val1, val2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("PIQP vs OSQP: larger LP (n=10)", {
  skip_if_not_installed("piqp")
  skip_if_not_installed("osqp")
  n <- 10
  set.seed(42)
  c_vec <- rnorm(n)
  A_mat <- matrix(rnorm(n * n), nrow = n)

  ## PIQP
  x1 <- Variable(n)
  prob1 <- Problem(Minimize(t(c_vec) %*% x1),
                   list(A_mat %*% x1 <= rep(10, n), x1 >= -5, x1 <= 5))
  val1 <- psolve(prob1, solver = "PIQP")
  expect_equal(status(prob1), "optimal")

  ## OSQP
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(t(c_vec) %*% x2),
                   list(A_mat %*% x2 <= rep(10, n), x2 >= -5, x2 <= 5))
  val2 <- psolve(prob2, solver = "OSQP")

  expect_equal(val1, val2, tolerance = 1e-2)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-2)
})
