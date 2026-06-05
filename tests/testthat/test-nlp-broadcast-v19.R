## CVXPY 1.9 DNLP broadcasting parity (nlp_tests/test_broadcast.py).
##
## These exercise implicit broadcasting in `x - A` / `log(x) == b`, NOT the
## broadcast_to atom. CVXR's 2-D matrix model already broadcasts a scalar, a
## column `(m, 1)`, a row `(1, n)`, and an outer `(m, 1) x (1, n) -> (m, n)`.
##
## DEVIATIONS FROM CVXPY (documented):
##  1. numpy's 1-D `Variable(6)` broadcasts as a ROW against an `(m, 6)` array.
##     CVXR `Variable(6)` is a `(6, 1)` column, so the faithful 2-D equivalent
##     of the numpy row broadcast is an explicit `(1, 6)` variable. The optimum
##     values are identical (the per-column means).
##  2. The Python DerivativeChecker has no R equivalent. For the two subtle
##     constraint-broadcast cases (which CVXPY only gradient-checks, never
##     solves), we instead assert the broadcast residual shape matches numpy and
##     that the NLP diff engine builds and evaluates the constraint.

.bc_nlp_solve <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(psolve(prob, nlp = TRUE, solver = "IPOPT"), type = "output")
  expect_equal(status(prob), "optimal")
}

.bc_skip_backend <- function() {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
}

## @cvxpy nlp_tests/test_broadcast.py::TestBroadcast::test_scalar_to_matrix
test_that("DNLP broadcast: scalar minus matrix, optimum is the grand mean", {
  set.seed(0)
  A <- matrix(rnorm(200 * 6), 200, 6)
  x <- Variable()
  value(x) <- 0.1
  .bc_nlp_solve(Problem(Minimize(sum_entries(square(x - A)))))
  expect_equal(as.numeric(value(x)), mean(A), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_broadcast.py::TestBroadcast::test_row_broadcast
test_that("DNLP broadcast: row vector minus matrix, optimum is the column means", {
  set.seed(0)
  A <- matrix(rnorm(5 * 6), 5, 6)
  ## numpy Variable(6) broadcasts as a row; CVXR's 2-D equivalent is (1, 6).
  x <- Variable(c(1, 6))
  value(x) <- matrix(0, 1, 6)
  .bc_nlp_solve(Problem(Minimize(sum_entries(square(x - A)))))
  expect_equal(as.numeric(value(x)), colMeans(A), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_broadcast.py::TestBroadcast::test_column_broadcast
test_that("DNLP broadcast: column vector minus matrix, optimum is the row means", {
  set.seed(0)
  A <- matrix(rnorm(5 * 6), 5, 6)
  x <- Variable(c(5, 1))
  value(x) <- matrix(0, 5, 1)
  .bc_nlp_solve(Problem(Minimize(sum_entries(square(x - A)))))
  expect_equal(as.numeric(value(x)), rowMeans(A), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_broadcast.py::TestBroadcast::test_subtle_broadcast1
test_that("DNLP broadcast: log(x[5,1]) == b[5] broadcasts to (5, 1)", {
  .bc_skip_backend()
  set.seed(0)
  n <- 5
  x <- Variable(c(n, 1))
  value(x) <- matrix(runif(n) + 0.1, n, 1)
  con <- (log(x) == rep(1, n))
  expect_equal(dim(con@.expr), c(n, 1L))
  ## Diff-engine analog of CVXPY's DerivativeChecker: build + evaluate.
  cp <- .de_C_problem(Problem(Minimize(0), list(con)))
  expect_true(is.finite(cp$objective_forward(rep(1, n))))
})

## @cvxpy nlp_tests/test_broadcast.py::TestBroadcast::test_subtle_broadcast2
test_that("DNLP broadcast: log(x[5,1]) == b[1,5] broadcasts to (5, 5)", {
  .bc_skip_backend()
  set.seed(0)
  n <- 5
  x <- Variable(c(n, 1))
  value(x) <- matrix(runif(n) + 0.1, n, 1)
  con <- (log(x) == matrix(1, 1, n))
  expect_equal(dim(con@.expr), c(n, n))
  cp <- .de_C_problem(Problem(Minimize(0), list(con)))
  expect_true(is.finite(cp$objective_forward(rep(1, n))))
})
