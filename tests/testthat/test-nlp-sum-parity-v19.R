## CVXPY 1.9 DNLP `sum` parity (nlp_tests/test_sum.py).
##
## sum is affine (smooth), so these problems are DNLP and solve through the NLP
## path. CVXPY runs them through IPOPT plus a Python DerivativeChecker; CVXR has
## no gradient checker, so these exercise the same models through
## psolve(..., nlp = TRUE) on the available R NLP backend. The two_sum tests use
## random data and (upstream) only verify derivatives, so here they verify the
## solve reaches an optimum. CVXPY axis 0 -> CVXR axis 2, axis 1 -> CVXR axis 1.

.sum_nlp_solve <- function(prob, ...) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(val <- psolve(prob, nlp = TRUE, ...), type = "output")
  expect_equal(status(prob), "optimal")
  invisible(val)
}

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_sum_without_axis
test_that("DNLP sum: minimize (sum(x) - 3)^2 with x <= 1", {
  x <- Variable(c(2, 1))
  value(x) <- matrix(0, 2, 1)
  prob <- Problem(Minimize((sum_entries(x) - 3)^2), list(x <= 1))
  .sum_nlp_solve(prob)
  expect_equal(value(x), matrix(1, 2, 1), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_sum_with_axis
test_that("DNLP sum: row-wise sum (axis = 1) least squares to 4", {
  X <- Variable(c(2, 3))
  value(X) <- matrix(0.5, 2, 3)
  prob <- Problem(Minimize(sum_entries((sum_entries(X, axis = 1) - 4)^2)),
                  list(X >= 0, X <= 1))
  .sum_nlp_solve(prob)
  expect_equal(value(X), matrix(1, 2, 3), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_sum_with_other_axis
test_that("DNLP sum: column-wise sum (axis = 2) least squares to 4", {
  X <- Variable(c(2, 3))
  value(X) <- matrix(0.5, 2, 3)
  prob <- Problem(Minimize(sum_entries((sum_entries(X, axis = 2) - 4)^2)),
                  list(X >= 0, X <= 1))
  .sum_nlp_solve(prob)
  expect_equal(value(X), matrix(1, 2, 3), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_two_sum_with_axis
test_that("DNLP sum: prod of row-wise sums of A %*% X solves", {
  set.seed(0)
  A <- matrix(runif(4 * 2), 4, 2)
  X <- Variable(c(2, 3))
  value(X) <- matrix(0.5, 2, 3)
  prob <- Problem(Minimize(prod_entries(sum_entries(A %*% X, axis = 1))),
                  list(X >= 0, X <= 1))
  .sum_nlp_solve(prob)
  expect_true(is.finite(value(prob)))
})

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_two_sum_with_other_axis
test_that("DNLP sum: prod of column-wise sums of A %*% X solves", {
  set.seed(0)
  A <- matrix(runif(4 * 2), 4, 2)
  X <- Variable(c(2, 3))
  value(X) <- matrix(0.5, 2, 3)
  prob <- Problem(Minimize(prod_entries(sum_entries(A %*% X, axis = 2))),
                  list(X >= 0, X <= 1))
  .sum_nlp_solve(prob)
  expect_true(is.finite(value(prob)))
})

## @cvxpy nlp_tests/test_sum.py::TestSumIPOPT::test_sum_matrix_arg
test_that("DNLP sum: minimize sum(A * T) over a box with nonneg A", {
  set.seed(0)
  n <- 40; m <- 20; k <- 4
  A <- matrix(runif(n * k), n, k) %*% matrix(runif(k * m), k, m)
  T <- Variable(c(n, m))
  value(T) <- matrix(1.5, n, m)
  prob <- Problem(Minimize(sum_entries(A * T)), list(T >= 1, T <= 2))
  .sum_nlp_solve(prob)
  expect_equal(value(T), matrix(1, n, m), tolerance = 1e-3)
})
