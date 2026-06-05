## CVXPY 1.9 DNLP scalar/matrix problem parity.
##
## CVXPY runs these through IPOPT plus a Python DerivativeChecker. CVXR's
## available R NLP backend is UNO/sparsediff, so these tests exercise the same
## scalar and matrix DNLP models through `psolve(..., nlp = TRUE)`.

nlp19_expect_optimal <- function(prob, expected = NULL, tol = 1e-3) {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  output <- capture.output(val <- psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
  expect_true(is.finite(val))
  if (!is.null(expected)) expect_equal(val, expected, tolerance = tol)
  invisible(val)
}

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_exp nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_logistic nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_logistic_matrix
test_that("DNLP scalar/matrix parity: exp and logistic problems solve", {
  x <- Variable()
  nlp19_expect_optimal(Problem(Minimize(exp(x)), list(x >= 4)),
                       exp(4), tol = 1e-4)

  y <- Variable()
  nlp19_expect_optimal(Problem(Minimize(logistic(y)), list(y >= 0.4)),
                       log1p(exp(0.4)), tol = 1e-4)

  Y <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Minimize(sum_entries(logistic(Y))), list(Y >= 0.4)),
                       6 * log1p(exp(0.4)), tol = 1e-4)

})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_entropy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_entropy_matrix
test_that("DNLP scalar/matrix parity: entropy problems solve", {
  x <- Variable()
  nlp19_expect_optimal(Problem(Maximize(entr(x)), list(x >= 0.1)),
                       1 / exp(1), tol = 1e-4)

  X <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Maximize(sum_entries(entr(X))), list(X >= 0.1)),
                       6 / exp(1), tol = 1e-4)
})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_KL nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_KL_matrix
test_that("DNLP scalar/matrix parity: KL problems solve", {
  p <- Variable()
  q <- Variable()
  nlp19_expect_optimal(
    Problem(Minimize(kl_div(p, q)), list(p >= 0.1, q >= 0.1, p + q == 2)),
    0, tol = 1e-5
  )
  expect_equal(as.numeric(value(p)), 1, tolerance = 1e-4)
  expect_equal(as.numeric(value(q)), 1, tolerance = 1e-4)

  X <- Variable(c(3, 3))
  Y <- Variable(c(3, 3))
  nlp19_expect_optimal(
    Problem(Minimize(sum_entries(kl_div(X, Y))),
            list(X >= 0.1, Y >= 0.1, sum_entries(X) + sum_entries(Y) == 6)),
    0, tol = 1e-5
  )
})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_power nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_power_matrix
test_that("DNLP scalar/matrix parity: power problems solve", {
  x <- Variable()
  nlp19_expect_optimal(Problem(Minimize(power(x, 3)), list(x >= 2)),
                       8, tol = 1e-4)

  X <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Minimize(sum_entries(power(X, 3))), list(X >= 2)),
                       6 * 8, tol = 1e-4)

  ## CVXPY's fractional-power tests also include minimizing `power(x, 0.6)`,
  ## which CVXR rejects as a concave minimization in its DNLP grammar. The
  ## convex fractional-power case is still exercised here as R-specific
  ## coverage, while the mixed upstream tests are classified N/A.
  y <- Variable()
  nlp19_expect_optimal(Problem(Minimize(power(y, 1.5)), list(y >= 2)),
                       2^1.5, tol = 1e-4)

  Y <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Minimize(sum_entries(power(Y, 1.5))), list(Y >= 2)),
                       6 * 2^1.5, tol = 1e-4)
})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_scalar_trig nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_matrix_trig nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_scalar_hyperbolic nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_matrix_hyperbolic
test_that("DNLP scalar/matrix parity: trig and hyperbolic problems solve", {
  for (atom in list(tan, sin, cos)) {
    x <- Variable()
    value(x) <- 1
    nlp19_expect_optimal(Problem(Minimize(atom(x)), list(x >= 0.1)))
  }

  for (atom in list(tan, sin, cos)) {
    X <- Variable(c(3, 2))
    value(X) <- matrix(1, 3, 2)
    nlp19_expect_optimal(Problem(Minimize(sum_entries(atom(X))), list(X >= 0.1)))
  }

  x <- Variable()
  value(x) <- 1
  nlp19_expect_optimal(Problem(Minimize(sinh(x)), list(x >= 0.1)))

  y <- Variable()
  value(y) <- 1
  nlp19_expect_optimal(Problem(Minimize(tanh(y)), list(y >= 0.1)))

  X <- Variable(c(3, 2))
  value(X) <- matrix(1, 3, 2)
  nlp19_expect_optimal(Problem(Minimize(sum_entries(sinh(X))), list(X >= 0.1)))

  Y <- Variable(c(3, 2))
  value(Y) <- matrix(1, 3, 2)
  nlp19_expect_optimal(Problem(Minimize(sum_entries(tanh(Y))), list(Y >= 0.1)))
})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_scalar_quad_form nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_scalar_quad_over_lin nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_matrix_quad_over_lin
test_that("DNLP scalar/matrix parity: quadratic forms solve", {
  x <- Variable()
  nlp19_expect_optimal(Problem(Minimize(quad_form(x, matrix(3, 1, 1))), list(x >= 1)),
                       3, tol = 1e-4)

  y <- Variable()
  z <- Variable()
  nlp19_expect_optimal(Problem(Minimize(quad_over_lin(y, z)), list(y >= 1, z <= 1)),
                       1, tol = 1e-4)

  X <- Variable(c(3, 2))
  s <- Variable()
  nlp19_expect_optimal(Problem(Minimize(quad_over_lin(X, s)), list(X >= 1, s <= 1)),
                       6, tol = 1e-4)
})

## @cvxpy nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_rel_entr_both_scalar_variables nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_rel_entr_matrix_variable_and_scalar_variable nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_rel_entr_scalar_variable_and_matrix_variable nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_rel_entr_both_matrix_variables nlp_tests/test_scalar_and_matrix_problems.py::TestScalarProblems::test_rel_entr_both_vector_variables
test_that("DNLP scalar/matrix parity: relative entropy broadcasting problems solve", {
  x <- Variable()
  y <- Variable()
  nlp19_expect_optimal(Problem(Minimize(rel_entr(x, y)),
                               list(x >= 0.1, y >= 0.1, x <= 2, y <= 2)))

  X <- Variable(c(3, 2))
  ys <- Variable()
  nlp19_expect_optimal(Problem(Minimize(sum_entries(rel_entr(X, ys))),
                               list(X >= 0.1, ys >= 0.1, X <= 2, ys <= 2)))

  xs <- Variable()
  Y <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Minimize(sum_entries(rel_entr(xs, Y))),
                               list(xs >= 0.1, Y >= 0.1, xs <= 2, Y <= 2)))

  A <- Variable(c(3, 2))
  B <- Variable(c(3, 2))
  nlp19_expect_optimal(Problem(Minimize(sum_entries(rel_entr(A, B))),
                               list(A >= 0.1, B >= 0.1, A <= 2, B <= 2)))

  u <- Variable(3)
  v <- Variable(3)
  nlp19_expect_optimal(Problem(Minimize(sum_entries(rel_entr(u, v))),
                               list(u >= 0.1, v >= 0.1, u <= 2, v <= 2)))
})
