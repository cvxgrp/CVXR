## CVXPY 1.9 NLP diff-engine stress tests for the normcdf atom
## (nlp_tests/stress_tests_diff_engine/test_normal_cdf.py).
##
## The R port runs the faithful DerivativeChecker mirror `nlp_derivative_check`
## then solves and asserts the optimum. cp.nlp.normcdf -> normcdf, x**2 -> x^2.

## @cvxpy nlp_tests/stress_tests_diff_engine/test_normal_cdf.py::TestNormalCdf::test_normcdf
test_that("DE normcdf: minimize normcdf((x - 1)^2) -> x = 1", {
  x <- Variable(); value(x) <- 0.5
  prob <- Problem(Minimize(normcdf((x - 1)^2)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)
  expect_equal(value(prob), 0.5, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_normal_cdf.py::TestNormalCdf::test_normcdf_with_quadratic
test_that("DE normcdf: maximize normcdf(x) - x^2 -> stationary", {
  lambda <- 1.0
  x <- Variable(); value(x) <- 0.2
  prob <- Problem(Maximize(normcdf(x) - lambda * x^2))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  xv <- as.numeric(value(x))
  grad <- exp(-0.5 * xv^2) / sqrt(2 * pi) - 2 * lambda * xv
  expect_equal(grad, 0.0, tolerance = 1e-4)
})
