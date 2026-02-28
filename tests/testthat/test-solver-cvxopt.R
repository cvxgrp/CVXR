## Tests for CVXOPT solver integration (backed by R cccp package)
## All tests are wrapped in skip_if_not_installed("cccp") so they
## pass cleanly when cccp is not available.

# ══════════════════════════════════════════════════════════════════
# Unit tests — constants and availability
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT_SOLVER constant is exported", {
  expect_equal(CVXR::CVXOPT_SOLVER, "CVXOPT")
})

## @cvxpy NONE
test_that("CVXOPT appears in installed_solvers when available", {
  skip_if_not_installed("cccp")
  solvers <- installed_solvers()
  expect_true("CVXOPT" %in% solvers)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT LP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT basic LP", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "CVXOPT")
  ## x = [2, 1], obj = 4
  expect_equal(val, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2.0, 1.0), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("CVXOPT LP with equality constraint", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] == 5, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, 5.0, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_1
test_that("CVXOPT lp_1 (CVXOPT classic example with duals)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_1()
  ## min -4x₀ - 5x₁, s.t. 2x₀ + x₁ ≤ 3, x₀ + 2x₁ ≤ 3, x ≥ 0
  ## Reference: x = [1, 1], obj = -9, duals = [1, 2, 0, 0]
  skip_if_not_installed("cccp")
  x <- Variable(2, name = "x")
  c1 <- (2 * x[1] + x[2] <= 3)
  c2 <- (x[1] + 2 * x[2] <= 3)
  c3 <- (x[1] >= 0)
  c4 <- (x[2] >= 0)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]), list(c1, c2, c3, c4))
  psolve(prob, solver = "CVXOPT")
  expect_equal(value(prob), -9.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), 0.0, tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(c4)), 0.0, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_2
test_that("CVXOPT lp_2 (bounds LP with duals)", {
  ## min x₀ + 0.5x₁, s.t. -100 ≤ x₀ ≤ -10, x₁ == 1
  ## Reference: x = [-100, 1], obj = -99.5, duals = [1, 0, -0.5]
  skip_if_not_installed("cccp")
  x <- Variable(2, name = "x")
  c1 <- (x[1] >= -100)
  c2 <- (x[1] <= -10)
  c3 <- (x[2] == 1)
  prob <- Problem(Minimize(x[1] + 0.5 * x[2]), list(c1, c2, c3))
  psolve(prob, solver = "CVXOPT")
  expect_equal(value(prob), -99.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(-100, 1), tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), -0.5, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CVXOPT Maximize LP", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, 10.0, tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT SOCP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_socp_0
test_that("CVXOPT socp_0 (norm + equality)", {
  ## min ||x||₂ + 1, s.t. x == 0
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2) + 1), list(x == 0))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("CVXOPT SOCP (norm minimization)", {
  skip_if_not_installed("cccp")
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 2)), list(sum_entries(x) == 1))
  val <- psolve(prob, solver = "CVXOPT")
  ## value = 1/sqrt(3), x = [1/3, 1/3, 1/3]
  expect_equal(val, 1 / sqrt(3), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_socp_1
test_that("CVXOPT socp_1 (SOC constraint with linear objective)", {
  ## min 3x₀ + 2x₁ + x₂, s.t. ||x||₂ ≤ y, x₀ + x₁ + 3x₂ ≥ 1, y ≤ 5
  ## Reference: obj ~ -13.55
  skip_if_not_installed("cccp")
  x <- Variable(3)
  y <- Variable()
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(p_norm(x, 2) <= y,
                       x[1] + x[2] + 3 * x[3] >= 1.0,
                       y <= 5))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, -13.5486, tolerance = 1e-1)
  expect_equal(as.numeric(value(y)), 5.0, tolerance = 1e-2)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT SDP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT SDP (min trace X s.t. X >> I)", {
  ## min trace(X) s.t. X - I is PSD, X is 2x2 symmetric
  ## Reference: X = I, value = 2.0, dual = I
  skip_if_not_installed("cccp")
  X <- Variable(c(2, 2), symmetric = TRUE)
  constr_psd <- (X %>>% diag(2))
  prob <- Problem(Minimize(matrix_trace(X)), list(constr_psd))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, 2.0, tolerance = 1e-3)
  Xv <- value(X)
  expect_equal(Xv[1,1], 1.0, tolerance = 1e-3)
  expect_equal(Xv[2,2], 1.0, tolerance = 1e-3)
  expect_equal(Xv[1,2], 0.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CVXOPT SDP with dual values", {
  ## min trace(X) s.t. X >> I, X 2x2 symmetric
  ## dual should be I (returned as n^2 flat vector, column-major)
  skip_if_not_installed("cccp")
  X <- Variable(c(2, 2), symmetric = TRUE)
  constr <- (X %>>% diag(2))
  prob <- Problem(Minimize(matrix_trace(X)), list(constr))
  psolve(prob, solver = "CVXOPT")
  dv <- dual_value(constr)
  ## Reshape flat dual to matrix
  dv_mat <- matrix(as.numeric(dv), 2, 2)
  expect_equal(dv_mat[1,1], 1.0, tolerance = 1e-2)
  expect_equal(dv_mat[2,2], 1.0, tolerance = 1e-2)
  expect_equal(dv_mat[1,2], 0.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("CVXOPT SDP 3x3", {
  ## min trace(X) s.t. X >> 2*I, X 3x3 symmetric
  ## Reference: X = 2I, value = 6.0
  skip_if_not_installed("cccp")
  X <- Variable(c(3, 3), symmetric = TRUE)
  constr <- (X %>>% (2 * diag(3)))
  prob <- Problem(Minimize(matrix_trace(X)), list(constr))
  val <- psolve(prob, solver = "CVXOPT")
  expect_equal(val, 6.0, tolerance = 1e-3)
  Xv <- value(X)
  expect_equal(diag(Xv), rep(2.0, 3), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT infeasible / unbounded
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT infeasible problem", {
  skip_if_not_installed("cccp")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x >= 5, x <= 1))
  psolve(prob, solver = "CVXOPT")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate",
                                    "solver_error"))
})

## @cvxpy NONE
test_that("CVXOPT unbounded problem", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 1))
  psolve(prob, solver = "CVXOPT")
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate",
                                    "solver_error"))
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT dual values
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT dual values (LP)", {
  skip_if_not_installed("cccp")
  ## min x + 2y  s.t.  x + y == 3, x >= 0, y >= 0
  ## Reference: value = 3, x = [3, 0]
  ## eq dual = -1, ineq duals = [~0, 1]
  x <- Variable(2)
  eq <- (x[1] + x[2] == 3)
  ineq <- (x >= 0)
  prob <- Problem(Minimize(x[1] + 2 * x[2]), list(eq, ineq))
  psolve(prob, solver = "CVXOPT")

  expect_equal(value(prob), 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(eq)), -1.0, tolerance = 1e-3)
  d <- as.numeric(dual_value(ineq))
  expect_equal(d[1], 0.0, tolerance = 1e-3)
  expect_equal(d[2], 1.0, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT solver options + standard params
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_options
test_that("CVXOPT solver_opts (tolerances)", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "CVXOPT",
                feastol = 1e-7, abstol = 1e-7, reltol = 1e-7,
                num_iter = 200L)
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("CVXOPT native solver_opts", {
  skip_if_not_installed("cccp")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "CVXOPT",
                maxiters = 50L, feastol = 1e-5)
  expect_equal(val, 1.0, tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity: CVXOPT vs Clarabel
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("CVXOPT vs Clarabel: LP value + primals + duals", {
  skip_if_not_installed("cccp")
  skip_if_not_installed("clarabel")

  ## CVXOPT (fresh objects)
  x1 <- Variable(2)
  eq1 <- (x1[1] + x1[2] == 3)
  ge1a <- (x1[1] >= 1)
  ge1b <- (x1[2] >= 0)
  prob1 <- Problem(Minimize(x1[1] + 2 * x1[2]), list(eq1, ge1a, ge1b))
  val1 <- psolve(prob1, solver = "CVXOPT")

  ## Clarabel (fresh objects)
  x2 <- Variable(2)
  eq2 <- (x2[1] + x2[2] == 3)
  ge2a <- (x2[1] >= 1)
  ge2b <- (x2[2] >= 0)
  prob2 <- Problem(Minimize(x2[1] + 2 * x2[2]), list(eq2, ge2a, ge2b))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(eq1)), as.numeric(dual_value(eq2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1a)), as.numeric(dual_value(ge2a)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1b)), as.numeric(dual_value(ge2b)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CVXOPT vs Clarabel: SOCP value + primals", {
  skip_if_not_installed("cccp")
  skip_if_not_installed("clarabel")

  ## CVXOPT
  x1 <- Variable(4)
  prob1 <- Problem(Minimize(p_norm(x1, 2)), list(sum_entries(x1) == 2))
  val1 <- psolve(prob1, solver = "CVXOPT")

  ## Clarabel
  x2 <- Variable(4)
  prob2 <- Problem(Minimize(p_norm(x2, 2)), list(sum_entries(x2) == 2))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CVXOPT vs Clarabel: SDP value + primals", {
  skip_if_not_installed("cccp")
  skip_if_not_installed("clarabel")

  ## CVXOPT
  X1 <- Variable(c(3, 3), symmetric = TRUE)
  prob1 <- Problem(Minimize(matrix_trace(X1)), list(X1 %>>% diag(3)))
  val1 <- psolve(prob1, solver = "CVXOPT")

  ## Clarabel
  X2 <- Variable(c(3, 3), symmetric = TRUE)
  prob2 <- Problem(Minimize(matrix_trace(X2)), list(X2 %>>% diag(3)))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-3)
  expect_equal(as.numeric(value(X1)), as.numeric(value(X2)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("CVXOPT vs Clarabel: larger LP (n=10)", {
  skip_if_not_installed("cccp")
  skip_if_not_installed("clarabel")
  n <- 10
  set.seed(42)
  c_vec <- rnorm(n)
  A_mat <- matrix(rnorm(n * n), nrow = n)

  ## CVXOPT
  x1 <- Variable(n)
  prob1 <- Problem(Minimize(t(c_vec) %*% x1),
                   list(A_mat %*% x1 <= rep(10, n), x1 >= -5, x1 <= 5))
  val1 <- psolve(prob1, solver = "CVXOPT")
  expect_equal(status(prob1), "optimal")

  ## Clarabel
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(t(c_vec) %*% x2),
                   list(A_mat %*% x2 <= rep(10, n), x2 >= -5, x2 <= 5))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-2)
})
