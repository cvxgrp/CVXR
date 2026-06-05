## CVXPY 1.9 DGP numeric variable bounds + infeasibility (test_dgp.py::TestDgp).
##
## The numeric-bounds counterpart of test_dgp_dpp.py (parametric bounds): a
## positive variable with numeric bounds solves through the DGP path with the
## bounds log-transformed into the log domain (Dgp2Dcp.variable_canon). The 3D
## axis cases (test_dgp_sum_3d_axis / test_dgp_pnorm_3d_axis) use ndim-3
## variables and are N/A for CVXR's 2D model (see notes/cvxpy_na_tests.txt).

## @cvxpy test_dgp.py::TestDgp::test_numeric_bounds
test_that("DGP numeric bounds: scalar two-sided", {
  x <- Variable(pos = TRUE, bounds = list(0.5, 5.0))
  psolve(Problem(Minimize(x)), gp = TRUE)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-4)

  x <- Variable(pos = TRUE, bounds = list(0.5, 5.0))
  psolve(Problem(Maximize(x)), gp = TRUE)
  expect_equal(as.numeric(value(x)), 5.0, tolerance = 1e-4)
})

## @cvxpy test_dgp.py::TestDgp::test_numeric_bounds_one_sided
test_that("DGP numeric bounds: one-sided", {
  x <- Variable(pos = TRUE, bounds = list(2.0, NULL))
  psolve(Problem(Minimize(x), list(x <= 10.0)), gp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)

  y <- Variable(pos = TRUE, bounds = list(NULL, 3.0))
  psolve(Problem(Maximize(y)), gp = TRUE)
  expect_equal(as.numeric(value(y)), 3.0, tolerance = 1e-4)
})

## @cvxpy test_dgp.py::TestDgp::test_numeric_bounds_vector
test_that("DGP numeric bounds: vector variable", {
  x <- Variable(3, pos = TRUE, bounds = list(0.5, 5.0))
  psolve(Problem(Minimize(sum_entries(x))), gp = TRUE)
  expect_equal(as.numeric(value(x)), rep(0.5, 3), tolerance = 1e-4)
})

## @cvxpy test_dgp.py::TestDgp::test_dgp_infeasible_no_crash
test_that("DGP infeasible problem reports infeasible without crashing", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x >= 2, x <= 1))
  psolve(prob, solver = "CLARABEL", gp = TRUE)
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})
