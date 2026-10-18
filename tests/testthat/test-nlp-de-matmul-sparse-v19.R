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

## TestSparseMatmulDispatch. The path is detected via the Lagrange Hessian
## lower-triangle nnz of sum_squares(A %*% x - b): a dense binding reports
## n(n+1)/2, the CSR binding the true A'A pattern. As in CVXPY 1.9.3 (D9 option
## 2, ADR D_193.5): a sparse Matrix factor always binds as CSR; a dense R matrix
## on the left binds as CSR below SPARSE_DENSITY_THRESHOLD (0.05) and dense at
## or above it; a dense matrix on the right always binds dense.
.de_mm_var <- function(n) {
  x <- Variable(n)
  value(x) <- rep(0, n)
  x
}

.de_hessian_lower_nnz <- function(prob, check = FALSE) {
  if (check) nlp_derivative_check(prob)
  skip_if_not_installed("sparsediff")
  ns <- asNamespace("CVXR")
  canon <- get("reduction_apply", envir = ns)(get("Dnlp2Smooth", envir = ns)(), prob)
  cp <- get(".de_C_problem", envir = ns)(canon[[1L]])
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  length(cp$hessian_sparsity()$rows)
}

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestSparseMatmulDispatch::test_below_threshold_routes_sparse
test_that("DE matmul-dispatch: a mostly-zero dense constant A gives the sparse Hessian", {
  set.seed(0); n <- 40
  full_lower <- n * (n + 1) / 2
  A_sp <- .de_sp_csr(Matrix::rsparsematrix(n, n, density = 0.05 / 5))
  A_dense <- as.matrix(A_sp)
  b <- rep(1, n)
  dense_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(A_dense %*% .de_mm_var(n) - b))), check = TRUE)
  csr_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(as_cvxr_expr(A_sp) %*% .de_mm_var(n) - b))))
  expect_lt(dense_nnz, full_lower)
  expect_equal(dense_nnz, csr_nnz)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestSparseMatmulDispatch::test_above_threshold_stays_dense
test_that("DE matmul-dispatch: a dense A above the threshold stays dense", {
  ## Block-diagonal, so the explicit CSR build is genuinely smaller and the
  ## full-triangle check is meaningful (density blk / n = 0.1 > 0.05).
  set.seed(0); n <- 40; blk <- 4
  expect_gt(blk / n, 0.05)
  full_lower <- n * (n + 1) / 2
  A_dense <- as.matrix(Matrix::bdiag(lapply(seq_len(n / blk),
                                            function(i) matrix(runif(blk^2), blk, blk))))
  b <- rep(1, n)
  dense_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(A_dense %*% .de_mm_var(n) - b))), check = TRUE)
  csr_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(as_cvxr_expr(.de_sp_csr(A_dense)) %*% .de_mm_var(n) - b))))
  expect_equal(dense_nnz, full_lower)
  expect_lt(csr_nnz, full_lower)
})

## CVXPY converters.py:58-65: a constant dense right factor always binds dense
## (no density switch on the right); a sparse one binds as CSR. No upstream test.
## @cvxpy NONE
test_that("DE matmul-dispatch: a dense right factor binds dense at any density", {
  set.seed(0); n <- 40
  full_lower <- n * (n + 1) / 2
  A_dense <- as.matrix(Matrix::rsparsematrix(n, n, density = 0.05 / 5))
  y <- function() { v <- Variable(c(1L, n)); value(v) <- matrix(0, 1L, n); v }
  b <- matrix(1, 1L, n)
  dense_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(y() %*% A_dense - b))), check = TRUE)
  csr_nnz <- .de_hessian_lower_nnz(
    Problem(Minimize(sum_squares(y() %*% as_cvxr_expr(.de_sp_csc(A_dense)) - b))))
  expect_equal(dense_nnz, full_lower)
  expect_lt(csr_nnz, full_lower)
})

## R-side saving (deliberate divergence from CVXPY's convert_expr,
## converters.py:112-119): a constant matmul factor or quad_form P that the
## converter reads by value() is not first converted as a child, so a sparse
## matrix is never densified (a 3000 x 3000 one was copied to 68.7 MB per solve).
## @cvxpy NONE
test_that("DE walker: a sparse constant factor or P is never densified", {
  skip_if_not_installed("sparsediff")
  set.seed(0); n <- 40
  A <- Matrix::rsparsematrix(n, n, density = 0.05) + Matrix::Diagonal(n)
  P <- Matrix::crossprod(A)
  densified <- character(0)
  orig <- get(".de_to_dense", envir = asNamespace("CVXR"))
  local_mocked_bindings(.de_to_dense = function(value) {
    if (inherits(value, "sparseMatrix")) densified <<- c(densified, class(value)[1L])
    orig(value)
  }, .package = "CVXR")
  x <- .de_mm_var(n)
  probs <- list(
    left  = Problem(Minimize(sum_squares(as_cvxr_expr(A) %*% x - 1))),
    right = Problem(Minimize(sum_squares(t(x) %*% as_cvxr_expr(A) - 1))),
    quad  = Problem(Minimize(quad_form(x, P))))
  for (nm in names(probs)) {
    densified <- character(0)
    .de_hessian_lower_nnz(probs[[nm]], check = TRUE)
    expect_identical(densified, character(0), info = nm)
  }
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_matmul_sparse.py::TestSparseMatmulDispatch::test_parametric_stays_dense
test_that("DE matmul-dispatch: a parametric A stays dense (full Hessian lower triangle)", {
  set.seed(0); n <- 40
  full_lower <- n * (n + 1) / 2
  A_dense <- as.matrix(Matrix::rsparsematrix(n, n, density = 0.05 / 5))
  P <- Parameter(c(n, n))
  value(P) <- A_dense
  prob <- Problem(Minimize(sum_squares(P %*% .de_mm_var(n) - rep(1, n))))
  expect_equal(.de_hessian_lower_nnz(prob), full_lower)
})
