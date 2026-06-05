## CVXPY 1.9 DGP variable bounds, numeric and parametric
## (test_dgp_dpp.py::TestDgpVariableBounds).
##
## A positive (DGP) variable may carry bounds that are numeric, Parameters, or
## log-log-affine Expressions. Dgp2Dcp.variable_canon log-transforms the bounds
## into the log domain and the downstream CvxAttr2Constr lowers them to
## constraints; parametric bounds canonicalize through the DGP tree without
## eagerly evaluating log(param.value) (CVXPY issue #3004). is_dpp(context) gates
## DGP-DPP vs DCP-DPP. SOLVER = CLARABEL, gp = TRUE.

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_parametric_bounds
test_that("DGP bounds: parametric lower/upper", {
  lb <- Parameter(pos = TRUE); value(lb) <- 2.0
  ub <- Parameter(pos = TRUE); value(ub) <- 5.0
  x <- Variable(pos = TRUE, bounds = list(lb, ub))
  psolve(Problem(Minimize(x)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)
  x <- Variable(pos = TRUE, bounds = list(lb, ub))
  psolve(Problem(Maximize(x)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 5.0, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_one_sided_parametric_bounds
test_that("DGP bounds: one-sided parametric", {
  p <- Parameter(pos = TRUE); value(p) <- 0.5
  x <- Variable(pos = TRUE, bounds = list(p, NULL))
  psolve(Problem(Minimize(x), list(x <= 10.0)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-4)

  q <- Parameter(pos = TRUE); value(q) <- 3.0
  y <- Variable(pos = TRUE, bounds = list(NULL, q))
  psolve(Problem(Maximize(y)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(y)), 3.0, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_mixed_numeric_and_parametric
test_that("DGP bounds: mixed numeric and parametric", {
  ub <- Parameter(pos = TRUE); value(ub) <- 4.0
  x <- Variable(pos = TRUE, bounds = list(0.5, ub))
  psolve(Problem(Minimize(x)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_re_solve
test_that("DGP bounds: re-solve after changing a parameter", {
  lb <- Parameter(pos = TRUE); value(lb) <- 2.0
  ub <- Parameter(pos = TRUE); value(ub) <- 10.0
  x <- Variable(pos = TRUE, bounds = list(lb, ub))
  prob <- Problem(Minimize(x))
  psolve(prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)
  value(lb) <- 3.0
  psolve(prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_is_dpp
test_that("DGP bounds: is_dpp context (dgp vs dcp)", {
  p <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE, bounds = list(p, NULL))
  expect_true(is_dpp(x, "dgp"))
  expect_true(is_dpp(Problem(Minimize(x), list(x <= 10.0)), "dgp"))

  ## Product of params: DPP in DGP (monomial), not in DCP.
  p2 <- Parameter(pos = TRUE)
  y <- Variable(pos = TRUE, bounds = list(p * p2, NULL))
  expect_true(is_dpp(y, "dgp"))
  expect_false(is_dpp(y, "dcp"))
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_enforce_dpp
test_that("DGP bounds: enforce_dpp on a DPP problem", {
  lb <- Parameter(pos = TRUE); value(lb) <- 2.0
  x <- Variable(pos = TRUE, bounds = list(lb, NULL))
  prob <- Problem(Minimize(x), list(x <= 100.0))
  psolve(prob, solver = "CLARABEL", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)
  value(lb) <- 7.0
  psolve(prob, solver = "CLARABEL", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(as.numeric(value(x)), 7.0, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_get_problem_data_without_param_values
## @v19-pending: get_problem_data on a DPP problem with an UNSET parametric bound [partial-now]
test_that("DGP bounds: get_problem_data without param values", {
  ## CVXPY (issue #3004) returns solver data even when the bound's Parameter has
  ## no value, because BOUNDED_VARIABLES-capable solvers carry the bound as a
  ## separate lower_bounds array. CVXR has not ported BOUNDED_VARIABLES, so the
  ## parametric bound lowers to a constraint whose numeric data requires the
  ## value, and problem_data raises "unspecified parameter". Revisit once
  ## bounds-as-solver-args land. The DGP canonicalization itself no longer
  ## eagerly evaluates log(param.value) (the core #3004 fix).
  skip("get_problem_data without param values needs BOUNDED_VARIABLES (bounds as solver args)")
  p <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE, bounds = list(p, NULL))
  d <- problem_data(Problem(Minimize(x), list(x <= 10.0)), solver = "CLARABEL", gp = TRUE)
  expect_false(is.null(d$data))
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_vector_and_matrix
test_that("DGP bounds: vector and matrix variables", {
  lb <- Parameter(pos = TRUE); value(lb) <- 2.0
  x <- Variable(3, pos = TRUE, bounds = list(lb, NULL))
  psolve(Problem(Minimize(sum_entries(x)), list(x <= 100.0)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), rep(2.0, 3), tolerance = 1e-4)

  X <- Variable(c(2L, 3L), pos = TRUE, bounds = list(lb, NULL))
  psolve(Problem(Minimize(sum_entries(X)), list(X <= 100.0)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(X)), rep(2.0, 6), tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_equal_parametric_bounds
test_that("DGP bounds: equal parametric bounds", {
  p <- Parameter(pos = TRUE); value(p) <- 3.0
  x <- Variable(pos = TRUE, bounds = list(p, p))
  psolve(Problem(Minimize(x)), solver = "CLARABEL", gp = TRUE)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-4)
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_solve_without_param_value_raises
test_that("DGP bounds: solving without a parameter value raises", {
  p <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE, bounds = list(p, NULL))
  expect_error(psolve(Problem(Minimize(x), list(x <= 10.0)), solver = "CLARABEL", gp = TRUE),
               "unspecified parameter")
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_non_log_log_affine_bounds_rejected
test_that("DGP bounds: non-log-log-affine bounds are not DGP-DPP", {
  p <- Parameter(3, pos = TRUE)
  ## norm(p) is log-log-convex but not log-log-affine.
  x <- Variable(pos = TRUE, bounds = list(p_norm(p, 2), NULL))
  expect_false(is_dpp(x, "dgp"))
  expect_false(is_dpp(Problem(Minimize(sqrt(x))), "dgp"))

  ## sum(p) is log-log-convex but not log-log-affine.
  y <- Variable(pos = TRUE, bounds = list(sum_entries(p), NULL))
  expect_false(is_dpp(y, "dgp"))
  expect_false(is_dpp(Problem(Minimize(sqrt(y))), "dgp"))

  ## p[1] (indexing) is log-log-affine, so this is DGP-DPP.
  z <- Variable(pos = TRUE, bounds = list(p[1], NULL))
  expect_true(is_dpp(z, "dgp"))
  expect_true(is_dpp(Problem(Minimize(sqrt(z))), "dgp"))
})
