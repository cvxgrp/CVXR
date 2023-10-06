context("test-g05-conic_solvers-glpk")
TOL <- 1e-6

test_that("test glpk lp 0", {
  StandardTestLPs.test_lp_0(solver = "GLPK")
})

test_that("test glpk lp 1", {
  StandardTestLPs.test_lp_1(solver = "GLPK")
})

test_that("test glpk lp 2", {
  StandardTestLPs.test_lp_2(solver = "GLPK")
})

test_that("test glpk lp 3", {
  StandardTestLPs.test_lp_3(solver = "GLPK")
})

test_that("test glpk lp 4", {
  StandardTestLPs.test_lp_4(solver = "GLPK")
})

test_that("test glpk lp 5", {
  StandardTestLPs.test_lp_5(solver = "GLPK")
})

test_that("test glpk lp 6", {
  StandardTestLPs.test_lp_6(solver = "GLPK")
})

test_that("test glpk mi lp 0", {
  StandardTestLPs.test_mi_lp_0(solver = "GLPK_MI")
})

test_that("test glpk mi lp 1", {
  StandardTestLPs.test_mi_lp_1(solver = "GLPK_MI")
})

test_that("test glpk mi lp 2", {
  StandardTestLPs.test_mi_lp_2(solver = "GLPK_MI")
})

test_that("test glpk mi lp 3", {
  StandardTestLPs.test_mi_lp_3(solver = "GLPK_MI")
})

test_that("test glpk mi lp 4", {
  StandardTestLPs.test_mi_lp_4(solver = "GLPK_MI")
})

test_that("test glpk mi lp 5", {
  StandardTestLPs.test_mi_lp_5(solver = "GLPK_MI")
})

test_that("test glpk options", {
  sth <- lp_1()
  require(cvxopt)
  # TODO: assert "tm_lim" not in cvxopt.glpk.options
  result <- solve(sth, solver = "GLPK", tm_lim = 100)
  # TODO: assert "tm_lim" not in cvxopt.glpk.options
  verify_objective(sth, tolerance = 1e-4)
  check_primal_feasibility(sth, tolerance = 1e-4)
  check_complementarity(sth, tolerance = 1e-4)
  check_dual_domains(sth, tolerance = 1e-4)
})

test_that("test glpk mi options", {
  sth <- mi_lp_1()
  require(cvxopt)
  # TODO: assert "tm_lim" not in cvxopt.glpk.options
  result <- solve(sth, solver = "GLPK_MI", tm_lim = 100, verbose = TRUE)
  # TODO: assert "tm_lim" not in cvxopt.glpk.options
  verify_objective(sth, tolerance = 1e-4)
  verify_primal_values(sth, tolerance = 1e-4)
})
