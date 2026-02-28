## Tests for FiniteSet constraint and Valinvec2mixedint reduction
## CVXPY reference: cvxpy/tests/test_valinvec2mixedint.py

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_1
test_that("FiniteSet basic: minimize sum, set = {1,3,5}", {
  skip_if_not_installed("Rglpk")
  x <- Variable(4)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1, 1), tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_1
test_that("FiniteSet ineq_form: same problem with ineq_form = TRUE", {
  skip_if_not_installed("Rglpk")
  x <- Variable(4)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5), ineq_form = TRUE)))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 4.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_2
test_that("FiniteSet single element", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  prob <- Problem(Minimize(x), list(FiniteSet(x, 7)))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 7.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_3
test_that("FiniteSet non-integer values", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  prob <- Problem(Maximize(x), list(FiniteSet(x, c(-1.5, 0.5, 2.5)), x <= 2))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  ## x <= 2 excludes 2.5, so best is 0.5
  expect_equal(opt_val, 0.5, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_4
test_that("FiniteSet duplicates in vec", {
  skip_if_not_installed("Rglpk")
  x <- Variable()
  prob <- Problem(Minimize(x), list(FiniteSet(x, c(1, 1, 1, 2, 2, 3, 3))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 1.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_5
test_that("FiniteSet 2D variable", {
  skip_if_not_installed("Rglpk")
  x <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(0, 1))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 0.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_6
test_that("FiniteSet affine expression", {
  skip_if_not_installed("Rglpk")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)),
    list(FiniteSet(2 * x[1] + 1, c(1, 3, 5)), x[2] == 0))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  ## 2*x[0]+1=1 → x[0]=0, x[1]=0 → sum=0
  expect_equal(opt_val, 0.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_non_affine_exception
test_that("FiniteSet non-affine errors", {
  x <- Variable()
  expect_error(FiniteSet(square(x), c(1, 2)), "affine")
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_7
test_that("FiniteSet two-element set", {
  skip_if_not_installed("Rglpk")
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 2))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 3.0, tolerance = 1e-4)
})

## @cvxpy test_valinvec2mixedint.py::TestFiniteSet::test_independent_entries
test_that("FiniteSet independent entries", {
  skip_if_not_installed("Rglpk")
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] + x[2] - x[3]),
    list(FiniteSet(x, c(0, 5, 10))))
  opt_val <- psolve(prob, solver = "GLPK_MI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, -10.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("FiniteSet is_dcp/is_dqcp", {
  x <- Variable(3)
  fs <- FiniteSet(x, c(1, 2, 3))
  expect_true(is_dcp(fs))
  expect_true(is_dqcp(fs))
})

## @cvxpy NONE
test_that("FiniteSet residual", {
  x <- Variable()
  fs <- FiniteSet(x, c(1, 3, 5))
  value(x) <- 2.0
  ## Closest to 2 in {1,3,5} is 1 (dist=1) or 3 (dist=1), so residual = 1
  expect_equal(residual(fs), 1.0)
})

## @cvxpy NONE
test_that("FiniteSet expr_name", {
  x <- Variable()
  fs <- FiniteSet(x, c(1, 2))
  nm <- expr_name(fs)
  expect_true(grepl("FiniteSet", nm))
})

## @cvxpy NONE
test_that("FiniteSet with Gurobi solver", {
  skip_if_not(requireNamespace("gurobi", quietly = TRUE))
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(x)), list(FiniteSet(x, c(1, 3, 5))))
  opt_val <- psolve(prob, solver = "GUROBI")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 3.0, tolerance = 1e-4)
})
