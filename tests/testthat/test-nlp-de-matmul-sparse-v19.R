## CVXPY 1.9 NLP diff-engine stress tests for matmul across matrix formats
## (nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py).
##
## Dense / CSR / CSC constants must produce identical oracles (and solutions).
## The R port runs the faithful DerivativeChecker mirror `nlp_derivative_check`
## at an in-domain point and, for the format-equivalence test, solves each format
## and asserts the solutions agree. Sparse constants are wrapped with
## as_cvxr_expr(); R's RsparseMatrix is the CSR analogue, dgCMatrix the CSC one.

.de_sp_csc <- function(M) Matrix::Matrix(M, sparse = TRUE)
.de_sp_csr <- function(M) as(Matrix::Matrix(M, sparse = TRUE), "RsparseMatrix")

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestMatmulDifferentFormats::test_dense_sparse_sparse
test_that("DE matmul-format: c'x s.t. A x == b, dense/CSR/CSC agree", {
  set.seed(0); n <- 10
  A <- matrix(runif(n * n), n, n)
  cc <- matrix(runif(n), n, 1)
  x0 <- matrix(runif(n), n, 1)
  b <- A %*% x0
  obj <- function(x) Minimize(t(cc) %*% x)
  sols <- list()
  for (nm in c("dense", "csr", "csc")) {
    x <- Variable(c(n, 1L), nonneg = TRUE)
    Ax <- switch(nm, dense = A, csr = as_cvxr_expr(.de_sp_csr(A)), csc = as_cvxr_expr(.de_sp_csc(A)))
    value(x) <- matrix(10, n, 1)
    prob <- Problem(obj(x), list(Ax %*% x == b))
    nlp_derivative_check(prob)
    .de_nlp_solve(prob)
    sols[[nm]] <- as.numeric(value(x))
  }
  expect_equal(sols$csr, sols$dense, tolerance = 1e-4)
  expect_equal(sols$csc, sols$dense, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestMatmulDifferentFormats::test_dense_left_matmul
test_that("DE matmul-format: sum_squares(A @ X - B) left", {
  set.seed(0); m <- 4; n <- 4
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  X <- Variable(c(n, n), nonneg = TRUE)
  value(X) <- matrix(runif(n * n), n, n)
  nlp_derivative_check(Problem(Minimize(sum_squares(A %*% X - B))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestMatmulDifferentFormats::test_dense_right_matmul
test_that("DE matmul-format: sum_squares(X @ A - B) right", {
  set.seed(0); m <- 4; n <- 4
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  X <- Variable(c(n, n), nonneg = TRUE)
  value(X) <- matrix(runif(n * n), n, n)
  nlp_derivative_check(Problem(Minimize(sum_squares(X %*% A - B))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestMatmulDifferentFormats::test_sparse_and_dense_matmul
test_that("DE matmul-format: sum_squares(A @ X @ C - B), C sparse", {
  set.seed(0); m <- 4; n <- 4
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  Cm <- matrix(runif(m * n), m, n); Cm[runif(m * n) < 0.5] <- 0; C <- as_cvxr_expr(.de_sp_csc(Cm))
  X <- Variable(c(n, n), nonneg = TRUE)
  value(X) <- matrix(runif(n * n), n, n)
  nlp_derivative_check(Problem(Minimize(sum_squares(A %*% X %*% C - B))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestMatmulDifferentFormats::test_sparse_and_dense_matmul2
test_that("DE matmul-format: sum_squares(C @ X @ A - B), C sparse", {
  set.seed(0); m <- 4; n <- 3
  A <- matrix(runif(n * m), n, m); B <- matrix(runif(m * m), m, m)
  Cm <- matrix(runif(m * n), m, n); Cm[runif(m * n) < 0.5] <- 0; C <- as_cvxr_expr(.de_sp_csc(Cm))
  X <- Variable(c(n, n), nonneg = TRUE)
  value(X) <- matrix(runif(n * n), n, n)
  nlp_derivative_check(Problem(Minimize(sum_squares(C %*% X %*% A - B))))
})
