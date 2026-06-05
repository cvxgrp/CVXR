## CVXPY 1.9 DNLP parametric-problem parity (nlp_tests/test_nlp_parameters.py).
##
## Parameters inside an NLP: a parametrized matrix factor in a matmul routes
## through the dense param-aware matmul (sd_left/right_matmul_dense fed the
## parameter's diff-engine node), and parametric scalars flow through
## multiply/promote. Each test solves the problem with two hardcoded constant
## values and again with a Parameter mutated between solves; the parametric and
## hardcoded solutions must agree.
##
## DEVIATIONS FROM CVXPY (documented):
##  1. CVXPY asserts the parametric and hardcoded solutions are bit-identical
##     (norm == 0). The parametric path uses a different diff-engine node (a
##     dense param matmul / param multiply) than the hardcoded path (a sparse
##     constant matmul), so we assert agreement to a tight tolerance.
##  2. CVXPY clears x.value (= None) before each solve and relies on IPOPT's
##     default initialization. The R NLP path requires an initial point, so the
##     same starting point is set before each solve (hardcoded and parametric
##     alike, so they converge to the same optimum).
##  3. The Python DerivativeChecker has no R equivalent and is dropped.
## test_parameter_broadcast (needs broadcast_to, not implemented in CVXR) and
## test_parameter_times_sparse_variable (needs sparse-pattern variables, not in
## CVXR's variable surface) are classified N/A in notes/cvxpy_na_tests.txt.

.par_solve <- function(prob, x, start) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  value(x) <- start
  capture.output(psolve(prob, nlp = TRUE, solver = "IPOPT"), type = "output")
  expect_equal(status(prob), "optimal")
  as.numeric(value(x))
}

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_entropy_maximization
test_that("DNLP parameters: entropy max with parametric A @ x <= b", {
  set.seed(0)
  m <- 10; n <- 5
  A1 <- matrix(runif(m * n), m, n); A2 <- matrix(runif(m * n), m, n)
  b1 <- rep(1, m); b2 <- rep(0.8, m)
  x <- Variable(n, nonneg = TRUE); start <- rep(1 / n, n)

  h1 <- .par_solve(Problem(Maximize(sum_entries(entr(x))), list(A1 %*% x <= b1, sum_entries(x) == 1)), x, start)
  h2 <- .par_solve(Problem(Maximize(sum_entries(entr(x))), list(A2 %*% x <= b2, sum_entries(x) == 1)), x, start)

  A <- Parameter(c(m, n)); b <- Parameter(m)
  prob <- Problem(Maximize(sum_entries(entr(x))), list(A %*% x <= b, sum_entries(x) == 1))
  value(A) <- A1; value(b) <- b1; p1 <- .par_solve(prob, x, start)
  value(A) <- A2; value(b) <- b2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_log_sum_exp
test_that("DNLP parameters: log_sum_exp(A @ x + b) with parametric A, b", {
  set.seed(0)
  m <- 10; n <- 5
  A1 <- matrix(rnorm(m * n), m, n); A2 <- matrix(rnorm(m * n), m, n)
  b1 <- rnorm(m); b2 <- rnorm(m)
  x <- Variable(n, bounds = list(-1, 1)); start <- rep(0, n)

  h1 <- .par_solve(Problem(Minimize(log_sum_exp(A1 %*% x + b1))), x, start)
  h2 <- .par_solve(Problem(Minimize(log_sum_exp(A2 %*% x + b2))), x, start)

  A <- Parameter(c(m, n)); b <- Parameter(m)
  prob <- Problem(Minimize(log_sum_exp(A %*% x + b)))
  value(A) <- A1; value(b) <- b1; p1 <- .par_solve(prob, x, start)
  value(A) <- A2; value(b) <- b2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_right_matmul
test_that("DNLP parameters: right matmul X @ A with parametric A, B", {
  set.seed(0)
  m <- 5; n <- 5; p <- 20
  A1 <- matrix(runif(n * p), n, p); A2 <- matrix(runif(n * p), n, p)
  B1 <- matrix(runif(m * p), m, p); B2 <- matrix(runif(m * p), m, p)
  X <- Variable(c(m, n), nonneg = TRUE); start <- matrix(0.5, m, n)

  h1 <- .par_solve(Problem(Minimize(sum_squares(X %*% A1 - B1))), X, start)
  h2 <- .par_solve(Problem(Minimize(sum_squares(X %*% A2 - B2))), X, start)

  A <- Parameter(c(n, p)); B <- Parameter(c(m, p))
  prob <- Problem(Minimize(sum_squares(X %*% A - B)))
  value(A) <- A1; value(B) <- B1; p1 <- .par_solve(prob, X, start)
  value(A) <- A2; value(B) <- B2; p2 <- .par_solve(prob, X, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_shared_across_expressions
test_that("DNLP parameters: A shared across objective and constraint", {
  set.seed(0)
  m <- 20; n <- 5
  A1 <- matrix(runif(m * n), m, n); A2 <- matrix(runif(m * n), m, n)
  b1 <- runif(m); b2 <- runif(m)
  x <- Variable(n, nonneg = TRUE); start <- rep(0.1, n)

  h1 <- .par_solve(Problem(Minimize(sum_squares(A1 %*% x - b1)), list(sum_entries(A1 %*% x) == 1)), x, start)
  h2 <- .par_solve(Problem(Minimize(sum_squares(A2 %*% x - b2)), list(sum_entries(A2 %*% x) == 1)), x, start)

  A <- Parameter(c(m, n)); b <- Parameter(m)
  prob <- Problem(Minimize(sum_squares(A %*% x - b)), list(sum_entries(A %*% x) == 1))
  value(A) <- A1; value(b) <- b1; p1 <- .par_solve(prob, x, start)
  value(A) <- A2; value(b) <- b2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_promote
test_that("DNLP parameters: multiply(promote(a, m), x) with parametric scalar a", {
  set.seed(0)
  m <- 10
  a1 <- 2.0; a2 <- 0.5
  b1 <- runif(m); b2 <- runif(m)
  x <- Variable(m, bounds = list(0, 10)); start <- rep(1, m)

  h1 <- .par_solve(Problem(Minimize(sum_squares(multiply(cvxr_promote(a1, c(m, 1L)), x) - b1))), x, start)
  h2 <- .par_solve(Problem(Minimize(sum_squares(multiply(cvxr_promote(a2, c(m, 1L)), x) - b2))), x, start)

  a <- Parameter(); b <- Parameter(m)
  prob <- Problem(Minimize(sum_squares(multiply(cvxr_promote(a, c(m, 1L)), x) - b)))
  value(a) <- a1; value(b) <- b1; p1 <- .par_solve(prob, x, start)
  value(a) <- a2; value(b) <- b2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_scalar_times_log
test_that("DNLP parameters: -sum(a * log(x)) with parametric scalar a", {
  set.seed(0)
  n <- 5
  a1 <- 2.0; a2 <- 0.5
  x <- Variable(n, nonneg = TRUE); start <- rep(1 / n, n)
  cons <- function(v) list(sum_entries(v) == 1)

  h1 <- .par_solve(Problem(Minimize(-sum_entries(a1 * log(x))), cons(x)), x, start)
  h2 <- .par_solve(Problem(Minimize(-sum_entries(a2 * log(x))), cons(x)), x, start)

  a <- Parameter()
  prob <- Problem(Minimize(-sum_entries(a * log(x))), cons(x))
  value(a) <- a1; p1 <- .par_solve(prob, x, start)
  value(a) <- a2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})

## @cvxpy nlp_tests/test_nlp_parameters.py::TestNlpParameters::test_parameter_scalar_times_log_reversed
test_that("DNLP parameters: -sum(log(x) * a) with parametric scalar a (reversed)", {
  set.seed(0)
  n <- 5
  a1 <- 2.0; a2 <- 0.5
  x <- Variable(n, nonneg = TRUE); start <- rep(1 / n, n)
  cons <- function(v) list(sum_entries(v) == 1)

  h1 <- .par_solve(Problem(Minimize(-sum_entries(log(x) * a1)), cons(x)), x, start)
  h2 <- .par_solve(Problem(Minimize(-sum_entries(log(x) * a2)), cons(x)), x, start)

  a <- Parameter()
  prob <- Problem(Minimize(-sum_entries(log(x) * a)), cons(x))
  value(a) <- a1; p1 <- .par_solve(prob, x, start)
  value(a) <- a2; p2 <- .par_solve(prob, x, start)
  expect_equal(p1, h1, tolerance = 1e-5)
  expect_equal(p2, h2, tolerance = 1e-5)
})
