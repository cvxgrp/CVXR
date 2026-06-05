## CVXPY 1.9 NLP diff-engine stress test for multiply across matrix formats
## (nlp_tests/stress_tests_diff_engine/test_multiply_sparse.py).
##
## sum(multiply(A, x)) over a box [-2, 2] drives x to -2 where A > 0 and to +2
## where A < 0, identically for dense / CSR / CSC A. The R port runs the faithful
## DerivativeChecker mirror then solves each format and asserts the sign pattern.

.de_sp_csc <- function(M) Matrix::Matrix(M, sparse = TRUE)
.de_sp_csr <- function(M) as(Matrix::Matrix(M, sparse = TRUE), "RsparseMatrix")

## @cvxpy nlp_tests/stress_tests_diff_engine/test_multiply_sparse.py::TestMultiplyDifferentFormats::test_dense_sparse_sparse
test_that("DE multiply-format: sum(A * x) box optimum, dense/CSR/CSC", {
  set.seed(0); n <- 5
  A <- matrix(runif(n * n) - 0.5, n, n)
  for (nm in c("dense", "csr", "csc")) {
    x <- Variable(c(n, n), bounds = list(-2, 2))
    Ae <- switch(nm, dense = A, csr = as_cvxr_expr(.de_sp_csr(A)), csc = as_cvxr_expr(.de_sp_csc(A)))
    value(x) <- matrix(runif(n * n), n, n)
    prob <- Problem(Minimize(sum_entries(Ae * x)))
    nlp_derivative_check(prob)
    .de_nlp_solve(prob)
    xv <- matrix(as.numeric(value(x)), n, n)
    expect_true(all(abs(xv[A > 0] - (-2)) < 1e-4))
    expect_true(all(abs(xv[A < 0] - 2) < 1e-4))
  }
})
