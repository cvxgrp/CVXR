## CVXPY 1.9 NLP diff-engine stress tests for the permuted-dense (PD) path
## (nlp_tests/stress_tests_diff_engine/test_permuted_dense.py).
##
## PD originates at left/right matmul of a dense (or sparse) constant with a leaf
## vector variable; these tests stress propagation of the PD Jacobian/Hessian
## structure through multiply / index / transpose / reshape-broadcast. CVXPY
## solves then runs DerivativeChecker; the R port runs the faithful mirror
## `nlp_derivative_check` at an in-domain point. sin/cos/logistic accept all
## reals, so any feasible point validates the engine. Sparse constants are
## wrapped with as_cvxr_expr(); CVXPY 0-based fancy indices become 1-based.

.de_sp <- function(M) Matrix::Matrix(M, sparse = TRUE)

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_pd_pd
test_that("DE PD: sin(A @ x) * cos(B @ y), dense/dense", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  nlp_derivative_check(Problem(Minimize(sum_entries(sin(A %*% x) * cos(B %*% y)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_pd_sparse
test_that("DE PD: sin(A @ x) * cos(B @ y), dense/sparse", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n)
  Bm <- matrix(runif(m * n), m, n); Bm[runif(m * n) < 0.5] <- 0; B <- as_cvxr_expr(.de_sp(Bm))
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  nlp_derivative_check(Problem(Minimize(sum_entries(sin(A %*% x) * cos(B %*% y)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_sparse_pd
test_that("DE PD: sin(A @ x) * cos(B @ y), sparse/dense", {
  set.seed(0); n <- 5; m <- 6
  Am <- matrix(runif(m * n), m, n); Am[runif(m * n) < 0.5] <- 0; A <- as_cvxr_expr(.de_sp(Am))
  B <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  nlp_derivative_check(Problem(Minimize(sum_entries(sin(A %*% x) * cos(B %*% y)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_pd_plain_var
test_that("DE PD: sin(A @ x) * cos(y), PD and plain var", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(m, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(m) * 2 - 1
  nlp_derivative_check(Problem(Minimize(sum_entries(sin(A %*% x) * cos(y)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_plain_var_pd
test_that("DE PD: sin(y) * cos(A @ x), plain var and PD", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(m, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(m) * 2 - 1
  nlp_derivative_check(Problem(Minimize(sum_entries(sin(y) * cos(A %*% x)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_index_propagation
test_that("DE PD: index into PD result (left matmul)", {
  set.seed(0); n <- 5; m <- 8
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  idx_A <- c(1L, 3L, 5L, 2L, 4L, 1L, 8L); idx_B <- c(1L, 5L, 3L, 4L, 2L, 1L, 8L)
  obj <- sum_entries(sin((A %*% x)[idx_A]) * cos((B %*% y)[idx_B]))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_transpose_propagation
test_that("DE PD: transpose of PD result (left matmul)", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  x <- Variable(c(n, 1L), bounds = list(-1, 1)); y <- Variable(c(n, 1L), bounds = list(-1, 1))
  value(x) <- matrix(runif(n) * 2 - 1, n, 1); value(y) <- matrix(runif(n) * 2 - 1, n, 1)
  obj <- sum_entries(sin(t(A %*% x)) * cos(t(B %*% y)))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_broadcast_propagation
test_that("DE PD: broadcast reshape of PD results (left matmul)", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  obj <- sum_entries(
    reshape_expr(sin(A %*% x), c(m, 1L), order = "F") *
      reshape_expr(cos(B %*% y), c(1L, m), order = "F"))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_deep_composition
test_that("DE PD: deep composition sin(A @ cos(B @ logistic(C @ x)))", {
  set.seed(0); n <- 5; m <- 10
  A <- matrix(runif(m * n), m, n); C <- matrix(runif(m * n), m, n)
  Bm <- matrix(runif(n * m), n, m); Bm[runif(n * m) < 0.5] <- 0; B <- as_cvxr_expr(.de_sp(Bm))
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  obj <- sum_entries(
    sin(A %*% cos(B %*% logistic(C %*% x))) *
      cos(A %*% cos(B %*% logistic(C %*% y))))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_pd_pd_right
test_that("DE PD: sin(x @ A) * cos(y @ B), dense/dense right", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(n * m), n, m); B <- matrix(runif(n * m), n, m)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  obj <- sum_entries(sin(t(x) %*% A) * cos(t(y) %*% B))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_multiply_pd_sparse_right
test_that("DE PD: sin(x @ A) * cos(y @ B), dense/sparse right", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(n * m), n, m)
  Bm <- matrix(runif(n * m), n, m); Bm[runif(n * m) < 0.5] <- 0; B <- as_cvxr_expr(.de_sp(Bm))
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  obj <- sum_entries(sin(t(x) %*% A) * cos(t(y) %*% B))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_index_propagation_right
test_that("DE PD: index into PD result (right matmul)", {
  set.seed(0); n <- 5; m <- 8
  A <- matrix(runif(n * m), n, m); B <- matrix(runif(n * m), n, m)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  idx_A <- c(1L, 3L, 5L, 2L, 4L, 1L, 8L); idx_B <- c(1L, 5L, 3L, 4L, 2L, 1L, 8L)
  obj <- sum_entries(sin((t(x) %*% A)[idx_A]) * cos((t(y) %*% B)[idx_B]))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_transpose_propagation_right
test_that("DE PD: transpose of PD result (right matmul)", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(n * m), n, m); B <- matrix(runif(n * m), n, m)
  x <- Variable(c(1L, n), bounds = list(-1, 1)); y <- Variable(c(1L, n), bounds = list(-1, 1))
  value(x) <- matrix(runif(n) * 2 - 1, 1, n); value(y) <- matrix(runif(n) * 2 - 1, 1, n)
  obj <- sum_entries(sin(t(x %*% A)) * cos(t(y %*% B)))
  nlp_derivative_check(Problem(Minimize(obj)))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_permuted_dense.py::TestPermutedDense::test_pd_broadcast_propagation_right
test_that("DE PD: broadcast reshape of PD results (right matmul)", {
  set.seed(0); n <- 5; m <- 6
  A <- matrix(runif(n * m), n, m); B <- matrix(runif(n * m), n, m)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n) * 2 - 1; value(y) <- runif(n) * 2 - 1
  obj <- sum_entries(
    reshape_expr(sin(t(x) %*% A), c(m, 1L), order = "F") *
      reshape_expr(cos(t(y) %*% B), c(1L, m), order = "F"))
  nlp_derivative_check(Problem(Minimize(obj)))
})
