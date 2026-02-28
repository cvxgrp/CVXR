## Tests for Phase C0: Complex number infrastructure prerequisites
## These verify the foundational fixes for complex number support.

## @cvxpy NONE
test_that("intf_convert handles complex scalars", {
  val <- intf_convert(1+2i)
  expect_true(is.matrix(val))
  expect_true(is.complex(val))
  expect_equal(dim(val), c(1L, 1L))
  expect_equal(val[1,1], 1+2i)
})

## @cvxpy NONE
test_that("intf_convert handles complex vectors", {
  val <- intf_convert(c(1+0i, 0+2i, 3+4i))
  expect_true(is.matrix(val))
  expect_true(is.complex(val))
  expect_equal(dim(val), c(3L, 1L))
  expect_equal(val[2,1], 0+2i)
})

## @cvxpy NONE
test_that("intf_convert handles complex matrices", {
  m <- matrix(c(1+0i, 0+1i, 0-1i, 2+0i), 2, 2)
  val <- intf_convert(m)
  expect_true(is.matrix(val))
  expect_true(is.complex(val))
  expect_equal(dim(val), c(2L, 2L))
  expect_equal(val, m)
})

## @cvxpy NONE
test_that("intf_convert preserves real matrices (no regression)", {
  m <- matrix(1:4, 2, 2)
  val <- intf_convert(m)
  expect_true(is.matrix(val))
  expect_false(is.complex(val))
})

## @cvxpy NONE
test_that("intf_sign returns unknown for complex values", {
  sgn <- intf_sign(matrix(1+2i, 1, 1))
  expect_false(sgn$is_nonneg)
  expect_false(sgn$is_nonpos)
})

## @cvxpy NONE
test_that("intf_sign still works for real values (no regression)", {
  sgn <- intf_sign(matrix(c(1, 2, 3), 3, 1))
  expect_true(sgn$is_nonneg)
  expect_false(sgn$is_nonpos)

  sgn2 <- intf_sign(matrix(c(-1, -2), 2, 1))
  expect_false(sgn2$is_nonneg)
  expect_true(sgn2$is_nonpos)
})

## @cvxpy NONE
test_that("intf_is_hermitian detects complex Hermitian matrices", {
  ## Hermitian: A == conj(A^T)
  H <- matrix(c(1+0i, 2+3i, 2-3i, 4+0i), 2, 2)
  res <- intf_is_hermitian(H)
  expect_true(res$is_hermitian)
  expect_false(res$is_symmetric)  # NOT complex-symmetric (A != A^T)
})

## @cvxpy NONE
test_that("intf_is_hermitian detects real symmetric (no regression)", {
  S <- matrix(c(1, 2, 2, 3), 2, 2)
  res <- intf_is_hermitian(S)
  expect_true(res$is_symmetric)
  expect_true(res$is_hermitian)
})

## @cvxpy NONE
test_that("intf_is_hermitian rejects non-Hermitian complex matrices", {
  ## Non-Hermitian complex matrix
  M <- matrix(c(1+0i, 2+3i, 4+5i, 6+0i), 2, 2)
  res <- intf_is_hermitian(M)
  expect_false(res$is_hermitian)
})

## @cvxpy NONE
test_that("intf_is_hermitian: complex-symmetric is NOT symmetric (CVXPY convention)", {
  ## CVXPY convention: "symmetric" always means "real symmetric".
  ## Complex matrices NEVER have is_symmetric=TRUE, even if A == A^T.
  ## CVXPY SOURCE: interface/matrix_utilities.py lines 278-284
  CS <- matrix(c(1+0i, 2+3i, 2+3i, 4+0i), 2, 2)
  res <- intf_is_hermitian(CS)
  expect_false(res$is_symmetric)   # complex → always FALSE
  expect_false(res$is_hermitian)   # not Hermitian (conj(A^T) != A)
})

## @cvxpy NONE
test_that("intf_is_skew_symmetric detects real skew-symmetric matrices", {
  S <- matrix(c(0, 2, -2, 0), 2, 2)
  expect_true(intf_is_skew_symmetric(S))
})

## @cvxpy NONE
test_that("intf_is_skew_symmetric rejects symmetric matrices", {
  S <- matrix(c(1, 2, 2, 3), 2, 2)
  expect_false(intf_is_skew_symmetric(S))
})

## @cvxpy NONE
test_that("intf_is_skew_symmetric returns FALSE for complex", {
  M <- matrix(c(0+0i, 1+1i, -1-1i, 0+0i), 2, 2)
  expect_false(intf_is_skew_symmetric(M))
})

# ── Constant complex detection ────────────────────────────────────────

## @cvxpy NONE
test_that("Constant(complex_scalar) is detected as complex", {
  c1 <- Constant(1+2i)
  expect_true(is_complex(c1))
  expect_false(is_imag(c1))
  expect_true(is_real(c1) == FALSE)  # is_real = !is_complex
})

## @cvxpy NONE
test_that("Constant(purely_imaginary) is detected as imaginary", {
  c1 <- Constant(0+3i)
  expect_true(is_complex(c1))
  expect_true(is_imag(c1))
})

## @cvxpy NONE
test_that("Constant(real) is not complex (no regression)", {
  c1 <- Constant(5)
  expect_false(is_complex(c1))
  expect_false(is_imag(c1))
  expect_true(is_real(c1))
})

## @cvxpy NONE
test_that("Constant complex matrix sign is unknown", {
  c1 <- Constant(matrix(c(1+2i, 3+4i), 2, 1))
  expect_false(is_nonneg(c1))
  expect_false(is_nonpos(c1))
})

## @cvxpy NONE
test_that("Constant(complex) is_pos returns FALSE", {
  c1 <- Constant(1+2i)
  expect_false(is_pos(c1))
})

# ── Constant Hermitian detection ──────────────────────────────────────

## @cvxpy NONE
test_that("Constant(Hermitian_matrix) is detected as Hermitian", {
  H <- matrix(c(1+0i, 2+3i, 2-3i, 4+0i), 2, 2)
  c1 <- Constant(H)
  expect_true(is_hermitian(c1))
  expect_true(is_complex(c1))
})

## @cvxpy NONE
test_that("Constant(real_symmetric) is detected as Hermitian", {
  S <- matrix(c(1, 2, 2, 3), 2, 2)
  c1 <- Constant(S)
  expect_true(is_hermitian(c1))
  expect_true(is_symmetric(c1))
})

## @cvxpy NONE
test_that("Constant(real_scalar) is Hermitian", {
  c1 <- Constant(5)
  expect_true(is_hermitian(c1))
})

## @cvxpy NONE
test_that("Constant(complex_scalar) is NOT Hermitian (complex scalar not real)", {
  c1 <- Constant(1+2i)
  ## Complex scalar: is_hermitian checks scalar && is_real, which is FALSE
  expect_false(is_hermitian(c1))
})

## @cvxpy NONE
test_that("Constant(non_Hermitian_complex) is not Hermitian", {
  M <- matrix(c(1+0i, 2+3i, 4+5i, 6+0i), 2, 2)
  c1 <- Constant(M)
  expect_false(is_hermitian(c1))
})

# ── Constant skew-symmetric detection ─────────────────────────────────

## @cvxpy NONE
test_that("Constant(skew_symmetric) detected", {
  S <- matrix(c(0, 3, -3, 0), 2, 2)
  c1 <- Constant(S)
  expect_true(is_skew_symmetric(c1))
})

## @cvxpy NONE
test_that("Constant(symmetric) not skew-symmetric", {
  S <- matrix(c(1, 2, 2, 3), 2, 2)
  c1 <- Constant(S)
  expect_false(is_skew_symmetric(c1))
})

## @cvxpy NONE
test_that("Constant(non-square) errors for skew-symmetric", {
  c1 <- Constant(matrix(1:6, 2, 3))
  expect_error(is_skew_symmetric(c1), "square")
})

## @cvxpy NONE
test_that("Constant(scalar) is skew-symmetric iff zero", {
  ## Scalar 0 is skew-symmetric: 0 + 0 == 0
  c0 <- Constant(0)
  expect_true(is_skew_symmetric(c0))
  ## Non-zero scalar is NOT skew-symmetric
  c1 <- Constant(5)
  expect_false(is_skew_symmetric(c1))
})

# ── Leaf.is_hermitian ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable(hermitian=TRUE) is_hermitian returns TRUE", {
  x <- Variable(c(3L, 3L), hermitian = TRUE)
  expect_true(is_hermitian(x))
  expect_true(is_complex(x))
})

## @cvxpy NONE
test_that("Variable(symmetric=TRUE) is_hermitian returns TRUE (real symmetric = Hermitian)", {
  x <- Variable(c(3L, 3L), symmetric = TRUE)
  expect_true(is_hermitian(x))
  expect_false(is_complex(x))
})

## @cvxpy NONE
test_that("Variable(PSD=TRUE) is_hermitian returns TRUE", {
  x <- Variable(c(3L, 3L), PSD = TRUE)
  expect_true(is_hermitian(x))
})

## @cvxpy NONE
test_that("Variable(NSD=TRUE) is_hermitian returns TRUE", {
  x <- Variable(c(3L, 3L), NSD = TRUE)
  expect_true(is_hermitian(x))
})

## @cvxpy NONE
test_that("plain Variable is not Hermitian (non-square)", {
  x <- Variable(c(3L, 1L))
  expect_false(is_hermitian(x))
})

## @cvxpy NONE
test_that("plain scalar Variable is Hermitian (real scalar)", {
  x <- Variable(c(1L, 1L))
  expect_true(is_hermitian(x))
})

# ── Expression.is_hermitian default ───────────────────────────────────

## @cvxpy NONE
test_that("Expression.is_hermitian requires is_real", {
  ## A real scalar expression should be Hermitian
  x <- Variable(c(1L, 1L))
  expect_true(is_hermitian(x))  # scalar, real → TRUE
})

# ── Expression.is_skew_symmetric default ──────────────────────────────

## @cvxpy NONE
test_that("Expression.is_skew_symmetric defaults to FALSE", {
  x <- Variable(c(3L, 3L))
  expect_false(is_skew_symmetric(x))
})
