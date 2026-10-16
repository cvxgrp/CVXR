## Tests for Gurobi solver integration
## All tests are wrapped in skip_if_not_installed("gurobi") so they
## pass cleanly when gurobi is not available.

# ══════════════════════════════════════════════════════════════════
# Unit tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GUROBI_SOLVER constant is exported", {
  expect_equal(CVXR::GUROBI_SOLVER, "GUROBI")
})

## @cvxpy NONE
test_that("Gurobi appears in installed_solvers when available", {
  skip_if_not_installed("gurobi")
  expect_true("GUROBI" %in% installed_solvers())
})

## @cvxpy NONE
test_that("Gurobi LP", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(3.0, 0.0), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi LP with equality constraint", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] == 5, x[1] >= 1, x[2] >= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 5.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi Maximize", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 10.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi QP", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 0.875, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.25, 0.75), ncol = 1), tolerance = 1e-3)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_var_bounds
test_that("Gurobi QP receives native variable bounds", {
  skip_if_not_installed("gurobi")
  x <- Variable(2, bounds = list(c(-1, 2), c(3, 4)))
  prob <- Problem(Minimize(sum_squares(x) + x[1] - 6 * x[2]))

  data <- problem_data(prob, solver = GUROBI_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], c(-1, 2))
  expect_equal(data[[UPPER_BOUNDS]], c(3, 4))
  expect_equal(data[[SD_DIMS]]@nonneg, 0L)

  val <- psolve(prob, solver = GUROBI_SOLVER)
  expect_equal(val, -9.25, tolerance = 1e-5)
  expect_equal(value(x), matrix(c(-0.5, 3), ncol = 1), tolerance = 1e-4)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_var_bounds
test_that("Gurobi QP updates parametric native bounds", {
  skip_if_not_installed("gurobi")
  lo <- Parameter(value = -1)
  hi <- Parameter(value = 2)
  x <- Variable(bounds = list(lo, hi))
  prob <- Problem(Minimize(square(x) + x))

  data <- problem_data(prob, solver = GUROBI_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]], -1)
  expect_equal(data[[UPPER_BOUNDS]], 2)

  val <- psolve(prob, solver = GUROBI_SOLVER)
  expect_equal(val, -0.25, tolerance = 1e-5)
  expect_equal(value(x), matrix(-0.5, ncol = 1), tolerance = 1e-4)

  value(lo) <- 0
  val2 <- psolve(prob, solver = GUROBI_SOLVER)
  expect_equal(val2, 0, tolerance = 1e-5)
  expect_equal(value(x), matrix(0, ncol = 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi SOCP via norm2", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, -sqrt(2), tolerance = 1e-4)
  expect_equal(value(x), matrix(c(-1/sqrt(2), -1/sqrt(2)), ncol = 1),
               tolerance = 1e-3)
})

## @cvxpy test_param_quad_prog.py::TestParamQuadProg::test_var_bounds
test_that("Gurobi conic path receives native variable bounds", {
  skip_if_not_installed("gurobi")
  x <- Variable(3, bounds = list(1, 3))
  prob <- Problem(Minimize(sum_entries(x)), list(cvxr_norm(x, 2) <= 10))

  data <- problem_data(prob, solver = GUROBI_SOLVER)$data
  expect_equal(data[[LOWER_BOUNDS]][seq_len(3)], rep(1, 3))
  expect_equal(data[[UPPER_BOUNDS]][seq_len(3)], rep(3, 3))
  expect_equal(tail(data[[LOWER_BOUNDS]], 1), -Inf)
  expect_equal(tail(data[[UPPER_BOUNDS]], 1), Inf)
  expect_true(length(data[[SD_DIMS]]@soc) > 0L)

  val <- psolve(prob, solver = GUROBI_SOLVER)
  expect_equal(val, 3, tolerance = 1e-5)
  expect_equal(value(x), matrix(rep(1, 3), ncol = 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi SOCP with multiple SOC constraints", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1]),
                  list(cvxr_norm(x, 2) <= 2,
                       x[2] >= 1))
  val <- psolve(prob, solver = "GUROBI")
  ## x2 = 1 (binding), x1 = -sqrt(4-1) = -sqrt(3) (on SOC boundary)
  expect_equal(val, -sqrt(3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi SOCP larger problem", {
  skip_if_not_installed("gurobi")
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] + x[2] + x[3]),
                  list(cvxr_norm(x, 2) <= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, -sqrt(3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi MIP LP (boolean)", {
  skip_if_not_installed("gurobi")
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(-x[1] - 2 * x[2]),
                  list(x[1] + x[2] <= 1.5))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, -2.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0, 1), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi MIP LP (integer)", {
  skip_if_not_installed("gurobi")
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(-x[1] - 2 * x[2]),
                  list(x[1] + x[2] <= 3, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, -6.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0, 3), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi MIQP (integer + quadratic)", {
  skip_if_not_installed("gurobi")
  ## Gurobi supports MIQP (unlike HiGHS)
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 2, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "GUROBI")
  ## Integer constraint: (x1,x2) in {(0,2),(1,1),(2,0)}
  ## obj at (0,2): 0+4+0=4, (1,1): 1+1+1=3, (2,0): 4+0+2=6
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(1, 1), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi MISOCP (boolean + SOC)", {
  skip_if_not_installed("gurobi")
  x <- Variable(2, boolean = TRUE)
  y <- Variable(2)
  prob <- Problem(Minimize(sum_entries(y)),
                  list(cvxr_norm(y, 2) <= 1,
                       y[1] >= x[1] - 1,
                       y[2] >= x[2] - 1,
                       x[1] + x[2] >= 1))
  val <- psolve(prob, solver = "GUROBI")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy NONE
test_that("Gurobi infeasible problem", {
  skip_if_not_installed("gurobi")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1, x <= -1))
  val <- psolve(prob, solver = "GUROBI")
  expect_true(status(prob) %in%
    c("infeasible", "infeasible_inaccurate", "infeasible_or_unbounded"))
  ## Value is +Inf for definite infeasible, NA for infeasible_or_unbounded
  expect_true(is.na(val) || (is.infinite(val) && val > 0))
})

## @cvxpy NONE
test_that("Gurobi infeasible with reoptimize", {
  skip_if_not_installed("gurobi")
  ## With reoptimize = TRUE, should get definitive answer
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1, x <= -1))
  val <- psolve(prob, solver = "GUROBI", reoptimize = TRUE)
  expect_true(status(prob) %in%
    c("infeasible", "infeasible_inaccurate"))
  expect_true(is.infinite(val) && val > 0)
})

## @cvxpy NONE
test_that("Gurobi unbounded problem", {
  skip_if_not_installed("gurobi")
  x <- Variable()
  prob <- Problem(Minimize(x))
  val <- psolve(prob, solver = "GUROBI")
  expect_true(status(prob) %in%
    c("unbounded", "unbounded_inaccurate", "infeasible_or_unbounded"))
  ## Value is -Inf for definite unbounded, NA for infeasible_or_unbounded
  expect_true(is.na(val) || (is.infinite(val) && val < 0))
})

## @cvxpy NONE
test_that("Gurobi solver options pass through", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x[1] + x[2] >= 4, x[1] >= 0, x[2] >= 0))
  ## Pass solver options (TimeLimit is harmless)
  val <- psolve(prob, solver = "GUROBI", TimeLimit = 100)
  expect_equal(val, 4.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi verbose output works", {
  skip_if_not_installed("gurobi")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 5))
  ## Should not error with verbose = TRUE
  val <- psolve(prob, solver = "GUROBI", verbose = TRUE)
  expect_equal(val, 5.0, tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# Dual variable tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Gurobi LP duals", {
  skip_if_not_installed("gurobi")
  ## min x1 + 2*x2  s.t.  x1+x2 >= 3, x1>=0, x2>=0
  ## At optimum x=(3,0): first constraint binding (dual nonzero),
  ## x1>=0 not binding (dual=0), x2>=0 binding (dual nonzero)
  x <- Variable(2)
  constr <- list(
    c1 <- (x[1] + x[2] >= 3),
    c2 <- (x[1] >= 0),
    c3 <- (x[2] >= 0)
  )
  prob <- Problem(Minimize(x[1] + 2 * x[2]), constr)
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 3.0, tolerance = 1e-4)
  ## Verify dual values exist and are nonneg for ineq constraints
  d1 <- dual_value(c1)
  expect_true(!is.null(d1))
  ## Binding constraint dual should be positive
  expect_true(d1 >= -1e-6)
})

## @cvxpy NONE
test_that("Gurobi QP duals", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  c1 <- (x[1] + x[2] == 1)
  c2 <- (x[1] >= 0)
  c3 <- (x[2] >= 0)
  prob <- Problem(Minimize(sum_squares(x) + x[1]), list(c1, c2, c3))
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 0.875, tolerance = 1e-4)
  ## Equality constraint should have a dual
  d1 <- dual_value(c1)
  expect_true(!is.null(d1))
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver parity tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Gurobi LP matches Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val_gurobi <- psolve(prob, solver = "GUROBI")

  ## Fresh problem for Clarabel (avoid cache leakage)
  x2 <- Variable(2)
  prob2 <- Problem(Minimize(x2[1] + 2 * x2[2]),
                   list(x2[1] + x2[2] >= 3, x2[1] >= 0, x2[2] >= 0))
  val_clarabel <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_gurobi, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi QP matches Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val_gurobi <- psolve(prob, solver = "GUROBI")

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_squares(x2) + x2[1]),
                   list(x2[1] + x2[2] == 1, x2[1] >= 0, x2[2] >= 0))
  val_clarabel <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_gurobi, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi SOCP matches Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val_gurobi <- psolve(prob, solver = "GUROBI")

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(x2[1] + x2[2]),
                   list(cvxr_norm(x2, 2) <= 1))
  val_clarabel <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_gurobi, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi Maximize matches Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Maximize(x[1] + x[2]),
                  list(x[1] + x[2] <= 10, x[1] >= 0, x[2] >= 0))
  val_gurobi <- psolve(prob, solver = "GUROBI")

  x2 <- Variable(2)
  prob2 <- Problem(Maximize(x2[1] + x2[2]),
                   list(x2[1] + x2[2] <= 10, x2[1] >= 0, x2[2] >= 0))
  val_clarabel <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_gurobi, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi MIP LP matches HiGHS", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("highs")
  x <- Variable(2, boolean = TRUE)
  prob <- Problem(Minimize(-x[1] - 2 * x[2]),
                  list(x[1] + x[2] <= 1.5))
  val_gurobi <- psolve(prob, solver = "GUROBI")

  x2 <- Variable(2, boolean = TRUE)
  prob2 <- Problem(Minimize(-x2[1] - 2 * x2[2]),
                   list(x2[1] + x2[2] <= 1.5))
  val_highs <- psolve(prob2, solver = "HIGHS")
  expect_equal(val_gurobi, val_highs, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Gurobi LP dual signs match Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")

  ## Fresh objects for each solver to avoid cache leakage
  ## Problem: min x1+2*x2 s.t. x1+x2=3, x1>=0, x2>=0
  x1 <- Variable(2)
  eq1 <- (x1[1] + x1[2] == 3)
  ge1 <- (x1[1] >= 0)
  ge2 <- (x1[2] >= 0)
  prob_g <- Problem(Minimize(x1[1] + 2 * x1[2]), list(eq1, ge1, ge2))
  psolve(prob_g, solver = "GUROBI")

  x2 <- Variable(2)
  eq2 <- (x2[1] + x2[2] == 3)
  ge3 <- (x2[1] >= 0)
  ge4 <- (x2[2] >= 0)
  prob_c <- Problem(Minimize(x2[1] + 2 * x2[2]), list(eq2, ge3, ge4))
  psolve(prob_c, solver = "CLARABEL")

  ## Equality constraint dual
  expect_equal(dual_value(eq1), dual_value(eq2), tolerance = 1e-3)
  ## Inequality constraint duals (sign should match)
  expect_equal(dual_value(ge1), dual_value(ge3), tolerance = 1e-3)
  expect_equal(dual_value(ge2), dual_value(ge4), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Gurobi QP dual signs match Clarabel", {
  skip_if_not_installed("gurobi")
  skip_if_not_installed("clarabel")

  x1 <- Variable(2)
  eq1 <- (x1[1] + x1[2] == 1)
  ge1 <- (x1[1] >= 0)
  ge2 <- (x1[2] >= 0)
  prob_g <- Problem(Minimize(sum_squares(x1) + x1[1]), list(eq1, ge1, ge2))
  psolve(prob_g, solver = "GUROBI")

  x2 <- Variable(2)
  eq2 <- (x2[1] + x2[2] == 1)
  ge3 <- (x2[1] >= 0)
  ge4 <- (x2[2] >= 0)
  prob_c <- Problem(Minimize(sum_squares(x2) + x2[1]), list(eq2, ge3, ge4))
  psolve(prob_c, solver = "CLARABEL")

  expect_equal(dual_value(eq1), dual_value(eq2), tolerance = 1e-2)
  expect_equal(dual_value(ge1), dual_value(ge3), tolerance = 1e-2)
  expect_equal(dual_value(ge2), dual_value(ge4), tolerance = 1e-2)
})

# ══════════════════════════════════════════════════════════════════
# CVXPY parity tests — standard test problems from solver_test_helpers
# Expected values verified via `uv run python` with CVXPY + gurobipy.
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_0
test_that("Gurobi parity: lp_0", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_0()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_1
test_that("Gurobi parity: lp_1 with duals", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_1()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
  ## Verify dual values
  cons <- constraints(sth$prob)
  for (i in seq_along(sth$con_duals)) {
    if (!is.null(sth$con_duals[[i]])) {
      expect_equal(as.numeric(dual_value(cons[[i]])),
                   sth$con_duals[[i]], tolerance = 1e-2)
    }
  }
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_2
test_that("Gurobi parity: lp_2 with duals", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_2()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
  cons <- constraints(sth$prob)
  for (i in seq_along(sth$con_duals)) {
    if (!is.null(sth$con_duals[[i]])) {
      expect_equal(as.numeric(dual_value(cons[[i]])),
                   sth$con_duals[[i]], tolerance = 1e-2)
    }
  }
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_3
test_that("Gurobi parity: lp_3 unbounded with reoptimize", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_3()
  val <- psolve(sth$prob, solver = "GUROBI", reoptimize = TRUE)
  expect_true(status(sth$prob) %in%
    c("unbounded", "unbounded_inaccurate", "infeasible_or_unbounded"))
  expect_true(is.na(val) || (is.infinite(val) && val < 0))
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_4
test_that("Gurobi parity: lp_4 infeasible with reoptimize", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_4()
  val <- psolve(sth$prob, solver = "GUROBI", reoptimize = TRUE)
  expect_true(status(sth$prob) %in%
    c("infeasible", "infeasible_inaccurate"))
  expect_true(is.infinite(val) && val > 0)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_5
test_that("Gurobi parity: lp_5", {
  skip_if_not_installed("gurobi")
  sth <- sth_lp_5()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("Gurobi parity: qp_0 with dual", {
  skip_if_not_installed("gurobi")
  sth <- sth_qp_0()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
  ## Dual for the binding constraint x >= 1
  cons <- constraints(sth$prob)
  expect_equal(as.numeric(dual_value(cons[[1]])),
               sth$con_duals[[1]], tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_0
test_that("Gurobi parity: socp_0", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_0()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_1
test_that("Gurobi parity: socp_1", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_1()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-3)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_2
test_that("Gurobi parity: socp_2", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_2()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-3)
})

## SOC duals: 1.9.3 plan defect 0.3. CVXR wrote the SOC coupling rows as
## [-A | -I][x; aux] = -b where CVXPY writes [A | I][x; aux] = b
## (gurobi_conif.py add_model_soc_constrs); the primal is the same, but
## Gurobi's Pi flips with the row, so every SOC dual came back negated.
## The tests above check objective/primal only, which is how it survived.
.gurobi_check_duals <- function(sth, tol) {
  for (i in seq_along(sth$con_duals)) {
    expected <- sth$con_duals[[i]]
    if (is.null(expected)) next
    actual <- dual_value(sth$prob@constraints[[i]])
    if (is.list(expected)) {
      for (k in seq_along(expected)) {
        expect_equal(as.numeric(actual[[k]]), as.numeric(expected[[k]]),
                     tolerance = tol, info = paste("constraint", i, "part", k))
      }
    } else {
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = tol,
                   info = paste("constraint", i))
    }
  }
}

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_4
test_that("Gurobi parity: socp_4 (three SOC blocks + equality) with duals", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_4()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-4)
  .gurobi_check_duals(sth, 1e-3)
})

## @cvxpy NONE
test_that("Gurobi SOC duals agree with the reference values on socp_1 and socp_2", {
  skip_if_not_installed("gurobi")
  for (make in list(sth_socp_1, sth_socp_2)) {
    sth <- make()
    psolve(sth$prob, solver = "GUROBI")
    .gurobi_check_duals(sth, 1e-3)
  }
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_3
test_that("Gurobi parity: socp_3 axis 0", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_3_ax0()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-3)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_3
test_that("Gurobi parity: socp_3 axis 1", {
  skip_if_not_installed("gurobi")
  sth <- sth_socp_3_ax1()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-3)
  x_var <- variables(sth$prob)[[1]]
  expect_equal(as.numeric(value(x_var)), sth$expect_x, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_lp_0
test_that("Gurobi parity: mi_lp_0", {
  skip_if_not_installed("gurobi")
  sth <- sth_mi_lp_0()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_lp_1
test_that("Gurobi parity: mi_lp_1", {
  skip_if_not_installed("gurobi")
  sth <- sth_mi_lp_1()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_lp_2
test_that("Gurobi parity: mi_lp_2 knapsack", {
  skip_if_not_installed("gurobi")
  sth <- sth_mi_lp_2()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_socp_1
test_that("Gurobi parity: mi_socp_1", {
  skip_if_not_installed("gurobi")
  sth <- sth_mi_socp_1()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_socp_2
test_that("Gurobi parity: mi_socp_2", {
  skip_if_not_installed("gurobi")
  sth <- sth_mi_socp_2()
  val <- psolve(sth$prob, solver = "GUROBI")
  expect_equal(val, sth$expect_obj, tolerance = 1e-4)
})

## Upstream reads each variable's VType back from the gurobipy model. The R gurobi
## result carries no model, so the types are read from the LP file written through
## save_file: its Binaries and Generals sections list the B and I columns.
## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mixed_integer_types
test_that("Gurobi: boolean, integer and continuous variables in one conic problem", {
  skip_if_not_installed("gurobi")
  z <- Variable(3, boolean = TRUE)
  w <- Variable(2, integer = TRUE)
  y <- Variable(2)
  objective <- Minimize(w[1] + 2 * w[2] - (3 * z[1] + 2 * z[2] + z[3]) + p_norm(y, 2))
  prob <- Problem(objective,
                  list(w >= -4, w <= 4, sum(z) <= 2, y >= 1, sum(w) >= -5))
  lp_file <- tempfile(fileext = ".lp")
  on.exit(unlink(lp_file), add = TRUE)
  val <- psolve(prob, solver = "GUROBI", save_file = lp_file)
  expect_equal(status(prob), "optimal")

  expect_equal(as.numeric(value(z)), c(1, 1, 0), tolerance = 1e-4)
  expect_equal(as.numeric(value(w)), c(-1, -4), tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), c(1, 1), tolerance = 1e-4)
  expect_equal(val, -14 + sqrt(2), tolerance = 1e-5)

  ## Each variable must reach Gurobi with the right type.
  lp <- trimws(readLines(lp_file))
  section <- function(name) {
    start <- match(name, lp)
    stop_at <- which(lp %in% c("Bounds", "Binaries", "Generals", "End"))
    stop_at <- min(stop_at[stop_at > start])
    unlist(strsplit(lp[seq.int(start + 1L, length.out = stop_at - start - 1L)], " +"))
  }
  expect_length(section("Binaries"), 3L)
  expect_length(section("Generals"), 2L)
})

## Gurobi minimizes x'Qx + c'x over the FULL Q. CVXR passed tril(P / 2), which
## counted every cross term P_ij x_i x_j once instead of twice, so a quad_form
## with an off-diagonal P solved a different problem and reported "optimal"
## (-0.612 here instead of -0.667). CVXPY passes the full 0.5 * P
## (gurobi_qpif.py:211-212). tril() also made a triangular dtCMatrix that the
## gurobi package then coerced back to dgCMatrix through a Matrix-deprecated
## route; the full matrix needs no coercion.
## @cvxpy NONE
test_that("Gurobi QP: off-diagonal quad_form matches the analytic optimum", {
  skip_if_not_installed("gurobi")
  old_opt <- options(Matrix.warnDeprecatedCoerce = 2L); on.exit(options(old_opt), add = TRUE)
  P <- matrix(c(2, 1.5, 1.5, 3), 2L)
  q <- c(-1, -2)
  x <- Variable(2)
  prob <- Problem(Minimize(0.5 * quad_form(x, P) + sum_entries(q * x)))
  ## Unconstrained minimizer x* = -P^{-1} q = (0, 2/3), optimum 0.5 q'x* = -2/3.
  x_star <- -solve(P, q)
  val <- psolve(prob, solver = "GUROBI")
  expect_equal(val, 0.5 * sum(q * x_star), tolerance = 1e-6)
  expect_equal(as.numeric(value(x)), x_star, tolerance = 1e-5)
  ## With x >= 0 binding: q2 = (1, -2) puts the unconstrained minimizer at
  ## (-1.6, 22/15), so the optimum is (0, 2/3) with dual (2, 0) on x >= 0.
  ## (With q, x1 >= 0 is active with a zero dual -- a degenerate point where
  ## interior-point x converges slowly; CVXPY + Gurobi gives the same x there.)
  q2 <- c(1, -2)
  x2 <- Variable(2)
  p2 <- Problem(Minimize(0.5 * quad_form(x2, P) + sum_entries(q2 * x2)), list(x2 >= 0))
  ref <- psolve(p2, solver = "CLARABEL")
  x3 <- Variable(2)
  p3 <- Problem(Minimize(0.5 * quad_form(x3, P) + sum_entries(q2 * x3)), list(x3 >= 0))
  expect_equal(psolve(p3, solver = "GUROBI"), ref, tolerance = 1e-6)
  expect_equal(as.numeric(value(x3)), c(0, 2 / 3), tolerance = 1e-5)
  expect_equal(as.numeric(dual_value(constraints(p3)[[1L]])), c(2, 0), tolerance = 1e-5)
})
