## CVXPY 1.9 NLP diff-engine stress tests for quad_form across matrix formats
## (nlp_tests/stress_tests_diff_engine/test_quad_form.py).
##
## A PSD P given as dense / CSR / CSC must yield the same QP solution, and a
## *parameterized* P must raise a clear error. The R port solves each format and
## asserts the solutions agree (plus the faithful DerivativeChecker mirror), and
## checks the parameterized-P error message.

.de_sp_csc <- function(M) Matrix::Matrix(M, sparse = TRUE)
.de_sp_csr <- function(M) as(Matrix::Matrix(M, sparse = TRUE), "RsparseMatrix")

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDifferentFormats::test_quad_form_dense_sparse_sparse
test_that("DE quad_form-format: (1/2) x'Px + q'x, column x, dense/CSR/CSC agree", {
  set.seed(1); m <- 15; n <- 10; p <- 5
  P <- matrix(rnorm(n * n), n, n); P <- t(P) %*% P
  q <- matrix(rnorm(n), n, 1)
  G <- matrix(rnorm(m * n), m, n); h <- G %*% matrix(rnorm(n), n, 1)
  A <- matrix(rnorm(p * n), p, n); b <- matrix(rnorm(p), p, 1)
  sols <- list()
  for (nm in c("dense", "csr", "csc")) {
    x <- Variable(c(n, 1L)); value(x) <- matrix(rnorm(n), n, 1)
    Pe <- switch(nm, dense = P, csr = as_cvxr_expr(.de_sp_csr(P)), csc = as_cvxr_expr(.de_sp_csc(P)))
    prob <- Problem(Minimize(0.5 * quad_form(x, Pe) + t(q) %*% x),
                    list(G %*% x <= h, A %*% x == b))
    nlp_derivative_check(prob)
    .de_nlp_solve(prob)
    sols[[nm]] <- as.numeric(value(x))
  }
  expect_equal(sols$csr, sols$dense, tolerance = 1e-4)
  expect_equal(sols$csc, sols$dense, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDifferentFormats::test_quad_form_dense_sparse_sparse_different_x
test_that("DE quad_form-format: (1/2) x'Px + q'x, 1-D x, dense/CSR/CSC agree", {
  set.seed(1); m <- 15; n <- 10; p <- 5
  P <- matrix(rnorm(n * n), n, n); P <- t(P) %*% P
  q <- matrix(rnorm(n), n, 1)
  G <- matrix(rnorm(m * n), m, n); h <- G %*% matrix(rnorm(n), n, 1)
  A <- matrix(rnorm(p * n), p, n); b <- matrix(rnorm(p), p, 1)
  sols <- list()
  for (nm in c("dense", "csr", "csc")) {
    x <- Variable(n); value(x) <- matrix(rnorm(n), n, 1)
    Pe <- switch(nm, dense = P, csr = as_cvxr_expr(.de_sp_csr(P)), csc = as_cvxr_expr(.de_sp_csc(P)))
    prob <- Problem(Minimize(0.5 * quad_form(x, Pe) + t(q) %*% x),
                    list(G %*% x <= h, A %*% x == b))
    nlp_derivative_check(prob)
    .de_nlp_solve(prob)
    sols[[nm]] <- as.numeric(value(x))
  }
  expect_equal(sols$csr, sols$dense, tolerance = 1e-4)
  expect_equal(sols$csc, sols$dense, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDifferentFormats::test_parameterized_quad_form_errors_clearly
test_that("DE quad_form-format: parameterized P errors clearly", {
  if (!.de_have_backend()) skip("No R NLP backend (sparsediff + ipopt/Uno)")
  set.seed(0); n <- 4
  x <- Variable(n, bounds = list(-1, 1))
  Pval <- matrix(runif(n * n), n, n); Pval <- Pval + t(Pval)
  P <- Parameter(c(n, n), symmetric = TRUE); value(P) <- Pval
  prob <- Problem(Minimize(quad_form(x, P)))
  expect_error(psolve(prob, nlp = TRUE, print_level = 0L), "parameterized P")
})
