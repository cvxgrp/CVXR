## CVXPY 1.9 NLP diff-engine stress tests for quad_form across matrix formats
## (nlp_tests/stress_tests_diff_engine/test_quad_form.py).
##
## A PSD P given as dense / CSR / CSC must yield the same QP solution. The R port
## solves each format and asserts the solutions agree (plus the faithful
## DerivativeChecker mirror). A parametric P is supported from CVXPY 1.9.3 on.

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

## TestQuadFormDiffEngine: derivative-engine checks for quad_form x'Px across P
## formats and a leaf or composed x (no NLP solver needed). `.de_spd` / `.de_qf_var`
## mirror upstream's `_spd` / `_var` (R's RNG, so the values differ; the checks
## are value-free).
.de_spd <- function(n, seed) {
  set.seed(seed); M <- matrix(rnorm(n * n), n, n)
  M %*% t(M) + n * diag(n)
}
.de_qf_var <- function(n, seed) {
  x <- Variable(n, bounds = list(-1, 1))
  set.seed(seed); value(x) <- runif(n, -0.9, 0.9)
  x
}

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_dense_P
test_that("DE quad_form: dense constant P over a leaf variable passes the derivative check", {
  prob <- Problem(Minimize(quad_form(.de_qf_var(6, 10), .de_spd(6, 0))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_sparse_P
test_that("DE quad_form: sparse constant P over a leaf variable passes the derivative check", {
  P <- as_cvxr_expr(.de_sp_csr(.de_spd(6, 1)))
  prob <- Problem(Minimize(quad_form(.de_qf_var(6, 11), P)))
  nlp_derivative_check(prob)
})

## CVXPY 1.9.3 (#3365) removed TestQuadFormDifferentFormats::test_parameterized_quad_form_errors_clearly:
## a parametric P is now supported by the diff engine, not rejected.
## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_parametric_P
test_that("DE quad_form: parametric P over a leaf variable passes the derivative check", {
  n <- 6
  P <- Parameter(c(n, n), PSD = TRUE); value(P) <- .de_spd(n, 2)
  prob <- Problem(Minimize(quad_form(.de_qf_var(n, 12), P)))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_dense_P_composition
test_that("DE quad_form: dense constant P over a sliced x[1:k] passes the derivative check", {
  x <- .de_qf_var(6, 13)
  prob <- Problem(Minimize(quad_form(x[1:3], .de_spd(3, 3))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_sparse_P_composition
test_that("DE quad_form: sparse constant P over a sliced x[1:k] passes the derivative check", {
  x <- .de_qf_var(6, 14)
  P <- as_cvxr_expr(.de_sp_csr(.de_spd(3, 4)))
  prob <- Problem(Minimize(quad_form(x[1:3], P)))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_dense_P_nonlinear_composition
test_that("DE quad_form: dense constant P over sin(x * x) passes the derivative check", {
  x <- .de_qf_var(6, 15)
  prob <- Problem(Minimize(quad_form(sin(multiply(x, x)), .de_spd(6, 5))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDiffEngine::test_quad_form_diag_param
test_that("DE quad_form: P = diag(p) for a vector Parameter p passes the derivative check", {
  n <- 6
  p <- Parameter(n, nonneg = TRUE); value(p) <- as.numeric(1:n)
  prob <- Problem(Minimize(quad_form(.de_qf_var(n, 16), diag(p))))
  nlp_derivative_check(prob)
})

## TestQuadFormDenseSparseDispatch: as in CVXPY 1.9.3 (registry.py:136-144; D9
## option 2, ADR D_193.5), a dense constant P below SPARSE_DENSITY_THRESHOLD (0.05)
## binds as CSR and one at or above it binds dense; a sparse Matrix P is always CSR.
## The R analogue of "routed to sparse" asserts the result that routing buys -- a
## diagonal Hessian block (n entries in the lower triangle), not n x n.
.de_hess_nnz <- function(prob) {
  ns <- asNamespace("CVXR")
  canon <- get("reduction_apply", envir = ns)(get("Dnlp2Smooth", envir = ns)(), prob)
  cp <- get(".de_C_problem", envir = ns)(canon[[1L]])
  cp$init_jacobian_coo()
  cp$init_hessian_coo_lower_tri()
  length(cp$hessian_sparsity()$rows)
}

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDenseSparseDispatch::test_dense_diagonal_P_routed_to_sparse
test_that("DE quad_form: a dense diagonal constant P gives a diagonal Hessian block", {
  skip_if_not_installed("sparsediff")
  n <- 40
  prob <- Problem(Minimize(quad_form(.de_qf_var(n, 21), diag(as.numeric(1:n)))))
  nlp_derivative_check(prob)
  expect_equal(.de_hess_nnz(prob), n)
})

## The above-threshold side of registry.py:136-144, which upstream does not test:
## a dense block-diagonal P at density 0.1 binds dense (full lower triangle), the
## same P as a sparse Matrix binds as CSR.
## @cvxpy NONE
test_that("DE quad_form: a dense P above the threshold binds dense", {
  skip_if_not_installed("sparsediff")
  set.seed(23); n <- 40; blk <- 4
  blocks <- lapply(seq_len(n / blk), function(i) {
    M <- matrix(rnorm(blk^2), blk, blk); crossprod(M) + diag(blk)
  })
  P <- as.matrix(Matrix::bdiag(blocks))
  dense_prob <- Problem(Minimize(quad_form(.de_qf_var(n, 24), P)))
  nlp_derivative_check(dense_prob)
  expect_equal(.de_hess_nnz(dense_prob), n * (n + 1) / 2)
  sparse_prob <- Problem(Minimize(quad_form(.de_qf_var(n, 24), .de_sp_csc(P))))
  expect_equal(.de_hess_nnz(sparse_prob), (n / blk) * blk * (blk + 1) / 2)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_quad_form.py::TestQuadFormDenseSparseDispatch::test_dense_diagonal_P_derivative_check
test_that("DE quad_form: a dense diagonal constant P passes the derivative check", {
  prob <- Problem(Minimize(quad_form(.de_qf_var(6, 22), diag(as.numeric(1:6)))))
  nlp_derivative_check(prob)
})
