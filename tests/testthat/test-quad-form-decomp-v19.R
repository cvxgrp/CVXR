## CVXPY 1.9 decomp_quad parity (test_quad_form.py::TestDecompQuad).
##
## decomp_quad(P) returns (scale, M1, M2) with P = scale * (M1 M1^H - M2 M2^H).
## CVXPY 1.9 computes this via LDL (Bunch-Kaufman) for dense indefinite matrices
## and sparse Cholesky for sparse PSD/NSD; CVXR uses a (Hermitian) eigenvalue
## decomposition. The decomposition *identity* is method-agnostic, so we assert
## it directly and drop CVXPY's LDL-internal checks (permutation, 2x2 pivots),
## which probe scipy's factorization rather than the decomp_quad contract.

.decomp_quad <- get("decomp_quad", envir = asNamespace("CVXR"))

## Mirror CVXPY's TestDecompQuad._check: P == scale * (M1 M1^H - M2 M2^H).
## decomp_quad emits the user-facing "nonconvex quad_form" warning iff P is
## indefinite (mixed-sign eigenvalues); `indefinite = TRUE` asserts that warning
## (CVXPY only suppresses it because Python lacks assert-and-consume).
.decomp_check <- function(P, tol = 1e-8, indefinite = FALSE) {
  if (indefinite) {
    expect_warning(res <- .decomp_quad(P), "nonconvex")
  } else {
    res <- .decomp_quad(P)
  }
  scale <- res[[1L]]; M1 <- res[[2L]]; M2 <- res[[3L]]
  Pd <- if (inherits(P, "Matrix")) as.matrix(P) else as.matrix(P)
  z <- matrix(0, nrow(Pd), ncol(Pd))
  pos <- if (length(M1)) scale * (M1 %*% Conj(t(M1))) else z
  neg <- if (length(M2)) scale * (M2 %*% Conj(t(M2))) else z
  expect_lt(max(Mod(Pd - (pos - neg))), tol)
}

.sp <- function(M) Matrix::Matrix(M, sparse = TRUE)

## @cvxpy test_quad_form.py::TestDecompQuad::test_indefinite
test_that("decomp_quad: dense indefinite (symmetric A + A^T)", {
  set.seed(0)
  for (n in c(3L, 10L, 30L)) {
    A <- matrix(rnorm(n * n), n, n)
    .decomp_check(A + t(A), indefinite = TRUE)
  }
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_2x2_blocks
## CVXPY checks that a zero diagonal forces Bunch-Kaufman 2x2 pivots; that is
## an LDL-internal property. CVXR uses eigendecomposition, so we keep only the
## decomposition-identity check on the same (zero-diagonal indefinite) matrices.
test_that("decomp_quad: indefinite with zero diagonal", {
  set.seed(0)
  for (n in c(6L, 10L, 30L)) {
    A <- matrix(rnorm(n * n), n, n)
    P <- A + t(A); diag(P) <- 0
    .decomp_check(P, indefinite = TRUE)
  }
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_complex_hermitian
test_that("decomp_quad: complex Hermitian (PSD and indefinite)", {
  set.seed(0); n <- 5
  A <- matrix(rnorm(n * n) + 1i * rnorm(n * n), n, n)
  .decomp_check(A %*% Conj(t(A)))                     # PSD
  .decomp_check(A + Conj(t(A)), indefinite = TRUE)    # indefinite Hermitian
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_sparse
test_that("decomp_quad: sparse scaled identity and signed diagonal", {
  .decomp_check(.sp(diag(5) * 3.0))                          # PSD
  .decomp_check(.sp(diag(c(1, -1, 2, -2, 3))), indefinite = TRUE)
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_sparse_nsd
test_that("decomp_quad: sparse NSD tridiagonal (negative Laplacian)", {
  n <- 6
  P <- matrix(0, n, n); diag(P) <- -4
  for (i in seq_len(n - 1L)) { P[i, i + 1L] <- 1; P[i + 1L, i] <- 1 }
  .decomp_check(.sp(P))
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_sparse_psd_rank_deficient
test_that("decomp_quad: sparse PSD/NSD rank-deficient", {
  set.seed(0)
  for (nk in list(c(4L, 2L), c(10L, 7L), c(20L, 5L))) {
    n <- nk[1L]; k <- nk[2L]
    B <- matrix(rnorm(n * k), n, k)
    .decomp_check(.sp(B %*% t(B)))      # PSD rank-deficient
    .decomp_check(.sp(-(B %*% t(B))))   # NSD rank-deficient
  }
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_sparse_diagonal_with_zeros
test_that("decomp_quad: sparse diagonal PSD with zero entries", {
  .decomp_check(.sp(diag(c(3.0, 0.0, 5.0, 0.0))))
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_sparse_zero
test_that("decomp_quad: sparse all-zero matrix", {
  res <- .decomp_quad(.sp(matrix(0, 4, 4)))
  expect_equal(res[[1L]], 0.0)
  expect_equal(length(res[[2L]]), 0L)
  expect_equal(length(res[[3L]]), 0L)
})
