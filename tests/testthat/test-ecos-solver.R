## Tests for ECOS and ECOS_BB solver integration
## All tests are wrapped in skip_if_not_installed("ECOSolveR") so they
## pass cleanly when ECOSolveR is not available.

# ══════════════════════════════════════════════════════════════════
# Unit tests — constants and availability
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS_SOLVER constant is exported", {
  expect_equal(CVXR::ECOS_SOLVER, "ECOS")
})

## @cvxpy NONE
test_that("ECOS_BB_SOLVER constant is exported", {
  expect_equal(CVXR::ECOS_BB_SOLVER, "ECOS_BB")
})

## @cvxpy NONE
test_that("ECOS appears in installed_solvers when available", {
  skip_if_not_installed("ECOSolveR")
  solvers <- installed_solvers()
  expect_true("ECOS" %in% solvers)
  expect_true("ECOS_BB" %in% solvers)
})

# ══════════════════════════════════════════════════════════════════
# ECOS LP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS basic LP", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "ECOS")
  ## CVXPY reference: value ≈ 4.0, x = [2, 1]
  expect_equal(val, 4.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(2.0, 1.0), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("ECOS LP with equality constraint", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] == 5, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "ECOS")
  expect_equal(val, 5.0, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_0
test_that("ECOS lp_0 (norm1 + equality)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_0()
  ## min ||x||₁ + 1, s.t. x == 0
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  val <- psolve(prob, solver = "ECOS")
  expect_equal(val, 1.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_1
test_that("ECOS lp_1 (CVXOPT example with duals)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_1()
  ## min -4x₀ - 5x₁, s.t. 2x₀ + x₁ ≤ 3, x₀ + 2x₁ ≤ 3, x ≥ 0
  ## CVXPY reference: x = [1, 1], obj = -9, duals = [1, 2, 0, 0]
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2, name = "x")
  c1 <- (2 * x[1] + x[2] <= 3)
  c2 <- (x[1] + 2 * x[2] <= 3)
  c3 <- (x[1] >= 0)
  c4 <- (x[2] >= 0)
  prob <- Problem(Minimize(-4 * x[1] - 5 * x[2]), list(c1, c2, c3, c4))
  psolve(prob, solver = "ECOS")
  expect_equal(value(prob), -9.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c4)), 0.0, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_2
test_that("ECOS lp_2 (bounds LP with duals)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_2()
  ## min x₀ + 0.5x₁, s.t. -100 ≤ x₀ ≤ -10, x₁ == 1
  ## CVXPY reference: x = [-100, 1], obj = -99.5, duals = [1, 0, -0.5]
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2, name = "x")
  c1 <- (x[1] >= -100)
  c2 <- (x[1] <= -10)
  c3 <- (x[2] == 1)
  prob <- Problem(Minimize(x[1] + 0.5 * x[2]), list(c1, c2, c3))
  psolve(prob, solver = "ECOS")
  expect_equal(value(prob), -99.5, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(-100, 1), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c1)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2)), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c3)), -0.5, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_3
test_that("ECOS lp_3 (unbounded with 5 vars)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_3()
  ## min sum(x), s.t. x <= 1  → unbounded
  skip_if_not_installed("ECOSolveR")
  x <- Variable(5)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 1))
  psolve(prob, solver = "ECOS")
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_4
test_that("ECOS lp_4 (infeasible with 5 vars)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_4()
  ## min sum(x), s.t. x <= 0, x >= 1  → infeasible
  skip_if_not_installed("ECOSolveR")
  x <- Variable(5)
  prob <- Problem(Minimize(sum_entries(x)), list(x <= 0, x >= 1))
  psolve(prob, solver = "ECOS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_lp_5
test_that("ECOS lp_5 (equality constraints, cross-check with Clarabel)", {
  ## CVXPY SOURCE: solver_test_helpers.py lp_5()
  ## Verify ECOS matches Clarabel on an LP with equality constraints
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")
  ## Fixed problem: min c^T x, s.t. Ax = b, x >= 0
  A_mat <- matrix(c(1, 0, 1, 0, 2,
                     0, 1, 0, 1, 0,
                     1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
  b_vec <- c(4, 3, 5)
  c_vec <- c(1, 2, 3, 1, 1)

  x1 <- Variable(5)
  prob1 <- Problem(Minimize(t(c_vec) %*% x1),
                   list(x1 >= 0, A_mat %*% x1 == b_vec))
  val1 <- psolve(prob1, solver = "ECOS")

  x2 <- Variable(5)
  prob2 <- Problem(Minimize(t(c_vec) %*% x2),
                   list(x2 >= 0, A_mat %*% x2 == b_vec))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(status(prob1), "optimal")
  expect_equal(val1, val2, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS Maximize", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "ECOS")
  expect_equal(val, 10.0, tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# ECOS SOCP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_socp_0
test_that("ECOS socp_0 (norm + equality)", {
  ## CVXPY SOURCE: solver_test_helpers.py socp_0()
  ## min ||x||₂ + 1, s.t. x == 0
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2) + 1), list(x == 0))
  val <- psolve(prob, solver = "ECOS")
  expect_equal(val, 1.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("ECOS SOCP (norm minimization)", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 2)), list(sum_entries(x) == 1))
  val <- psolve(prob, solver = "ECOS")
  ## CVXPY reference: value ≈ 1/sqrt(3) ≈ 0.5774, x = [1/3, 1/3, 1/3]
  expect_equal(val, 1 / sqrt(3), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_socp_1
test_that("ECOS socp_1 (SOC constraint with linear objective)", {
  ## CVXPY SOURCE: solver_test_helpers.py socp_1()
  ## min 3x₀ + 2x₁ + x₂, s.t. ||x||₂ ≤ y, x₀ + x₁ + 3x₂ ≥ 1, y ≤ 5
  ## CVXPY reference: obj ≈ -13.5486
  skip_if_not_installed("ECOSolveR")
  x <- Variable(3)
  y <- Variable()
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(p_norm(x, 2) <= y,
                       x[1] + x[2] + 3 * x[3] >= 1.0,
                       y <= 5))
  val <- psolve(prob, solver = "ECOS")
  expect_equal(val, -13.5486, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 5.0, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# ECOS ExpCone tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS ExpCone (entropy)", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  prob <- Problem(Maximize(sum_entries(entr(x))),
                  list(sum_entries(x) == 1, x >= 0.01))
  val <- psolve(prob, solver = "ECOS")
  ## CVXPY reference: value ≈ ln(2) ≈ 0.6931, x = [0.5, 0.5]
  expect_equal(val, log(2), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_expcone_1
test_that("ECOS expcone_1 (direct ExpCone constraint)", {
  ## CVXPY SOURCE: solver_test_helpers.py expcone_1()
  ## min 3x₀ + 2x₁ + x₂, s.t. 0.1 ≤ Σx ≤ 1, x ≥ 0, ExpCone(x₂, x₁, x₀)
  ## Tests expcone_permutor (critical for ECOS dual permutation)
  ## CVXPY reference: obj ≈ 0.2353, x ≈ [0.0546, 0.0261, 0.0193]
  skip_if_not_installed("ECOSolveR")
  x <- Variable(c(3, 1))
  cone_con <- ExpCone(x[3], x[2], x[1])
  c1 <- (sum_entries(x) <= 1.0)
  c2 <- (sum_entries(x) >= 0.1)
  c3 <- (x >= 0)
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(c1, c2, c3, cone_con))
  psolve(prob, solver = "ECOS")

  expect_equal(value(prob), 0.23535, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 0.05463, tolerance = 1e-2)
  expect_equal(xv[2], 0.02609, tolerance = 1e-2)
  expect_equal(xv[3], 0.01928, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_exp_soc_1
test_that("ECOS exp_soc_1 (risk-parity portfolio: ExpCone + SOC)", {
  ## CVXPY SOURCE: solver_test_helpers.py expcone_socp_1()
  ## Mixed exponential cone + SOC (risk-parity portfolio)
  ## CVXPY reference: obj ≈ 4.0751
  skip_if_not_installed("ECOSolveR")
  sigma <- matrix(c(1.83, 1.79, 3.22,
                     1.79, 2.18, 3.18,
                     3.22, 3.18, 8.69), nrow = 3, byrow = TRUE)
  L <- t(chol(sigma))  # lower Cholesky
  c_val <- 0.75
  t_var <- Variable(name = "t")
  x <- Variable(3, name = "x")
  s <- Variable(3, name = "s")
  e <- rep(1, 3)

  con1 <- (p_norm(t(L) %*% x, 2) <= t_var)
  con2 <- ExpCone(s, Constant(matrix(e, ncol = 1)), x)
  prob <- Problem(Minimize(t_var - c_val * sum_entries(s)), list(con1, con2))
  psolve(prob, solver = "ECOS")

  expect_equal(value(prob), 4.0751, tolerance = 1e-2)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 0.576, tolerance = 1e-2)
  expect_equal(xv[2], 0.543, tolerance = 1e-2)
  expect_equal(xv[3], 0.280, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("ECOS ExpCone (log_sum_exp)", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(3)
  prob <- Problem(Minimize(log_sum_exp(x)), list(sum_entries(x) == 0))
  val <- psolve(prob, solver = "ECOS")
  ## log(3 * exp(0)) = log(3) when all x_i = 0 (by symmetry)
  expect_equal(val, log(3), tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# ECOS dual values
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS dual values", {
  skip_if_not_installed("ECOSolveR")
  ## min x + 2y  s.t.  x + y == 3, x >= 0, y >= 0
  ## CVXPY reference: value = 3, x = [3, 0]
  ## eq dual = -1, ineq duals = [~0, 1]
  x <- Variable(2)
  eq <- (x[1] + x[2] == 3)
  ineq <- (x >= 0)
  prob <- Problem(Minimize(x[1] + 2 * x[2]), list(eq, ineq))
  psolve(prob, solver = "ECOS")

  expect_equal(value(prob), 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(eq)), -1.0, tolerance = 1e-4)
  d <- as.numeric(dual_value(ineq))
  expect_equal(d[1], 0.0, tolerance = 1e-4)
  expect_equal(d[2], 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("ECOS duals match Clarabel", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")

  ## Solve with ECOS (fresh objects)
  x1 <- Variable(name = "x1")
  y1 <- Variable(name = "y1")
  eq1 <- (x1 + y1 == 3)
  ineq1a <- (x1 >= 1)
  ineq1b <- (y1 >= 0)
  prob1 <- Problem(Minimize(x1 + 2 * y1), list(eq1, ineq1a, ineq1b))
  psolve(prob1, solver = "ECOS")

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

# ══════════════════════════════════════════════════════════════════
# ECOS infeasible / unbounded
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS infeasible problem", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x >= 5, x <= 1))
  psolve(prob, solver = "ECOS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy NONE
test_that("ECOS unbounded problem", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(x), list(x <= 10))
  psolve(prob, solver = "ECOS")
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

# ══════════════════════════════════════════════════════════════════
# ECOS solver options
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_options
test_that("ECOS solver_opts (all tolerance params)", {
  ## CVXPY SOURCE: test_conic_solvers.py test_ecos_options()
  ## Tests all 7 ECOS control parameters + warm_start
  skip_if_not_installed("ECOSolveR")
  EPS <- 1e-4
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  ## Run twice — second call exercises warm_start
  for (i in seq_len(2)) {
    val <- psolve(prob, solver = "ECOS",
                  feastol = EPS, abstol = EPS, reltol = EPS,
                  feastol_inacc = EPS, abstol_inacc = EPS,
                  reltol_inacc = EPS, max_iters = 20L)
  }
  expect_equal(val, 1.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# ECOS_BB MIP tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS_BB boolean variable", {
  skip_if_not_installed("ECOSolveR")
  ## min sum(x) s.t. x[1] + x[2] >= 1.5, x boolean
  ## Need at least 2 of 3 true → opt = 2
  x <- Variable(3, boolean = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x[1] + x[2] >= 1.5))
  val <- psolve(prob, solver = "ECOS_BB")
  ## CVXPY reference: value ≈ 2.0, x = [1, 1, 0]
  expect_equal(val, 2.0, tolerance = 1e-3)
  xv <- round(as.numeric(value(x)))
  expect_equal(sum(xv), 2L)
})

## @cvxpy NONE
test_that("ECOS_BB integer variable", {
  skip_if_not_installed("ECOSolveR")
  ## min x + 2y s.t. x + y >= 3, x,y integer >= 0
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(x[1] + 2 * x[2]), list(sum_entries(x) >= 3, x >= 0))
  val <- psolve(prob, solver = "ECOS_BB")
  ## CVXPY reference: value ≈ 3.0, x = [3, 0]
  expect_equal(val, 3.0, tolerance = 1e-2)
  xv <- round(as.numeric(value(x)))
  expect_equal(xv[1] + 2L * xv[2], 3L)
})

## @cvxpy NONE
test_that("ECOS_BB no duals for MIP", {
  skip_if_not_installed("ECOSolveR")
  x <- Variable(integer = TRUE, name = "x_int")
  constr <- (x >= 1)
  prob <- Problem(Minimize(x), list(constr))
  psolve(prob, solver = "ECOS_BB")
  ## MIP problems should not have dual values
  expect_null(dual_value(constr))
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_mi_lp_0
test_that("ECOS_BB mi_lp_0 (boolean with equality)", {
  ## CVXPY SOURCE: solver_test_helpers.py mi_lp_0()
  ## min ||x||₁ + 1, s.t. x == bool_var, bool_var == 0
  skip_if_not_installed("ECOSolveR")
  x <- Variable(2)
  bool_var <- Variable(boolean = TRUE)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0),
                  list(x == bool_var, bool_var == 0))
  val <- psolve(prob, solver = "ECOS_BB")
  expect_equal(val, 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-3)
  expect_equal(round(as.numeric(value(bool_var))), 0)
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_mi_lp_3
test_that("ECOS_BB mi_lp_3 (infeasible boolean MIP)", {
  ## CVXPY SOURCE: solver_test_helpers.py mi_lp_3()
  ## 4 boolean vars, sum == 2, but pairwise constraints make it infeasible
  skip_if_not_installed("ECOSolveR")
  x <- Variable(4, boolean = TRUE)
  prob <- Problem(Maximize(Constant(1)),
                  list(x[1] + x[2] + x[3] + x[4] <= 2,
                       x[1] + x[2] + x[3] + x[4] >= 2,
                       x[1] + x[2] <= 1,
                       x[1] + x[3] <= 1,
                       x[1] + x[4] <= 1,
                       x[3] + x[4] <= 1,
                       x[2] + x[4] <= 1,
                       x[2] + x[3] <= 1))
  psolve(prob, solver = "ECOS_BB")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy NONE
test_that("ECOS_BB solver_opts (mi params)", {
  ## Test ECOS_BB-specific mixed-integer options
  skip_if_not_installed("ECOSolveR")
  x <- Variable(3, boolean = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x[1] + x[2] >= 1.5))
  val <- psolve(prob, solver = "ECOS_BB",
                mi_max_iters = 1000L,
                mi_abs_eps = 1e-6,
                mi_rel_eps = 1e-6)
  expect_equal(round(val), 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS_BB integer LP where integrality matters", {
  skip_if_not_installed("ECOSolveR")
  ## max x s.t. x <= 3.5, x integer, x >= 0
  ## Continuous optimum: x = 3.5, integer optimum: x = 3
  x <- Variable(integer = TRUE, name = "x_int")
  prob <- Problem(Maximize(x), list(x <= 3.5, x >= 0))
  val <- psolve(prob, solver = "ECOS_BB")
  expect_equal(val, 3.0, tolerance = 1e-3)
  expect_equal(round(as.numeric(value(x))), 3L)
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity: ECOS vs Clarabel
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS vs Clarabel: LP value + primals + duals", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")
  ## min x₀ + 2x₁ + 3x₂  s.t. x₀+x₁+x₂ == 6, x ≥ 1, x₀+x₁ ≤ 5

  ## ECOS (fresh objects)
  x1 <- Variable(3)
  eq1 <- (x1[1] + x1[2] + x1[3] == 6)
  ge1 <- (x1 >= 1)
  le1 <- (x1[1] + x1[2] <= 5)
  prob1 <- Problem(Minimize(x1[1] + 2 * x1[2] + 3 * x1[3]), list(eq1, ge1, le1))
  val1 <- psolve(prob1, solver = "ECOS")

  ## Clarabel (fresh objects)
  x2 <- Variable(3)
  eq2 <- (x2[1] + x2[2] + x2[3] == 6)
  ge2 <- (x2 >= 1)
  le2 <- (x2[1] + x2[2] <= 5)
  prob2 <- Problem(Minimize(x2[1] + 2 * x2[2] + 3 * x2[3]), list(eq2, ge2, le2))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-4)
  expect_equal(as.numeric(dual_value(eq1)), as.numeric(dual_value(eq2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1)), as.numeric(dual_value(ge2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(le1)), as.numeric(dual_value(le2)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS vs Clarabel: SOCP value + primals", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")

  ## ECOS
  x1 <- Variable(4)
  prob1 <- Problem(Minimize(p_norm(x1, 2)), list(sum_entries(x1) == 2))
  val1 <- psolve(prob1, solver = "ECOS")

  ## Clarabel
  x2 <- Variable(4)
  prob2 <- Problem(Minimize(p_norm(x2, 2)), list(sum_entries(x2) == 2))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("ECOS vs Clarabel: ExpCone value + primals", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")

  ## ECOS
  x1 <- Variable(3, pos = TRUE)
  prob1 <- Problem(Maximize(sum_entries(entr(x1))),
                   list(sum_entries(x1) == 1))
  val1 <- psolve(prob1, solver = "ECOS")

  ## Clarabel
  x2 <- Variable(3, pos = TRUE)
  prob2 <- Problem(Maximize(sum_entries(entr(x2))),
                   list(sum_entries(x2) == 1))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("ECOS vs Clarabel: log_sum_exp", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")

  ## ECOS
  x1 <- Variable(4)
  prob1 <- Problem(Minimize(log_sum_exp(x1)), list(sum_entries(x1) == 1))
  val1 <- psolve(prob1, solver = "ECOS")

  ## Clarabel
  x2 <- Variable(4)
  prob2 <- Problem(Minimize(log_sum_exp(x2)), list(sum_entries(x2) == 1))
  val2 <- psolve(prob2, solver = "CLARABEL")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity: ECOS vs HiGHS
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS vs HiGHS: LP value + primals + duals", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("highs")
  ## Non-degenerate LP with unique optimal vertex:
  ## min -4x₀ - 5x₁  s.t. 2x₀+x₁ ≤ 3, x₀+2x₁ ≤ 3, x ≥ 0
  ## Unique optimum: x = [1, 1], obj = -9

  ## ECOS (fresh objects)
  x1 <- Variable(2)
  c1a <- (2 * x1[1] + x1[2] <= 3)
  c2a <- (x1[1] + 2 * x1[2] <= 3)
  c3a <- (x1 >= 0)
  prob1 <- Problem(Minimize(-4 * x1[1] - 5 * x1[2]), list(c1a, c2a, c3a))
  val1 <- psolve(prob1, solver = "ECOS")

  ## HiGHS (fresh objects)
  x2 <- Variable(2)
  c1b <- (2 * x2[1] + x2[2] <= 3)
  c2b <- (x2[1] + 2 * x2[2] <= 3)
  c3b <- (x2 >= 0)
  prob2 <- Problem(Minimize(-4 * x2[1] - 5 * x2[2]), list(c1b, c2b, c3b))
  val2 <- psolve(prob2, solver = "HIGHS")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c1a)), as.numeric(dual_value(c1b)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(c2a)), as.numeric(dual_value(c2b)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS vs HiGHS: LP with equality + inequality (unique optimum)", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("highs")
  ## min x + 2y  s.t. x + y == 5, x >= 2, y >= 1
  ## Unique optimum: x = 4, y = 1 (y at lower bound), obj = 6

  ## ECOS
  x1 <- Variable(2)
  eq1 <- (x1[1] + x1[2] == 5)
  ge1a <- (x1[1] >= 2)
  ge1b <- (x1[2] >= 1)
  prob1 <- Problem(Minimize(x1[1] + 2 * x1[2]), list(eq1, ge1a, ge1b))
  val1 <- psolve(prob1, solver = "ECOS")

  ## HiGHS
  x2 <- Variable(2)
  eq2 <- (x2[1] + x2[2] == 5)
  ge2a <- (x2[1] >= 2)
  ge2b <- (x2[2] >= 1)
  prob2 <- Problem(Minimize(x2[1] + 2 * x2[2]), list(eq2, ge2a, ge2b))
  val2 <- psolve(prob2, solver = "HIGHS")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(eq1)), as.numeric(dual_value(eq2)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1a)), as.numeric(dual_value(ge2a)), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(ge1b)), as.numeric(dual_value(ge2b)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS vs HiGHS: Maximize LP", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("highs")
  ## max 5x + 4y  s.t. x + y <= 10, 2x + y <= 14, x,y >= 0

  ## ECOS
  x1 <- Variable(2)
  prob1 <- Problem(Maximize(5 * x1[1] + 4 * x1[2]),
                   list(x1[1] + x1[2] <= 10,
                        2 * x1[1] + x1[2] <= 14,
                        x1 >= 0))
  val1 <- psolve(prob1, solver = "ECOS")

  ## HiGHS
  x2 <- Variable(2)
  prob2 <- Problem(Maximize(5 * x2[1] + 4 * x2[2]),
                   list(x2[1] + x2[2] <= 10,
                        2 * x2[1] + x2[2] <= 14,
                        x2 >= 0))
  val2 <- psolve(prob2, solver = "HIGHS")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity: ECOS vs all three (Clarabel, HiGHS, SCS)
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ECOS vs Clarabel vs HiGHS: larger LP (n=20)", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")
  skip_if_not_installed("highs")
  n <- 20
  set.seed(42)
  c_vec <- rnorm(n)
  A_mat <- matrix(rnorm(n * n), nrow = n)

  ## ECOS
  x1 <- Variable(n)
  prob1 <- Problem(Minimize(t(c_vec) %*% x1),
                   list(A_mat %*% x1 <= rep(10, n), x1 >= -5, x1 <= 5))
  val1 <- psolve(prob1, solver = "ECOS")
  expect_equal(status(prob1), "optimal")

  ## Clarabel
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(t(c_vec) %*% x2),
                   list(A_mat %*% x2 <= rep(10, n), x2 >= -5, x2 <= 5))
  val2 <- psolve(prob2, solver = "CLARABEL")

  ## HiGHS
  x3 <- Variable(n)
  prob3 <- Problem(Minimize(t(c_vec) %*% x3),
                   list(A_mat %*% x3 <= rep(10, n), x3 >= -5, x3 <= 5))
  val3 <- psolve(prob3, solver = "HIGHS")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(val1, val3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x3)), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("ECOS vs Clarabel vs SCS: SOCP (n=8)", {
  skip_if_not_installed("ECOSolveR")
  skip_if_not_installed("clarabel")
  skip_if_not_installed("scs")
  set.seed(7)
  n <- 8
  c_vec <- rnorm(n)
  A_mat <- matrix(rnorm(4 * n), nrow = 4)
  b_vec <- rnorm(4)

  ## ECOS
  x1 <- Variable(n)
  prob1 <- Problem(Minimize(p_norm(x1, 2)),
                   list(A_mat %*% x1 == b_vec))
  val1 <- psolve(prob1, solver = "ECOS")

  ## Clarabel
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(p_norm(x2, 2)),
                   list(A_mat %*% x2 == b_vec))
  val2 <- psolve(prob2, solver = "CLARABEL")

  ## SCS
  x3 <- Variable(n)
  prob3 <- Problem(Minimize(p_norm(x3, 2)),
                   list(A_mat %*% x3 == b_vec))
  val3 <- psolve(prob3, solver = "SCS")

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(val1, val3, tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
})
