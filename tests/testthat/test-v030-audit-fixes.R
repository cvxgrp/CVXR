## Tests for v0.3.0 audit fixes
## Covers: DGP/log-log chain, class hierarchy, sign propagation, etc.

library(testthat)
library(CVXR)

## ═══════════════════════════════════════════════════════════════════
## DGP / Log-log Chain (B1-B4, H7-H13, H25)
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("is_log_log_constant: positive constants are log-log constant", {
  c1 <- Constant(5)
  expect_true(is_log_log_constant(c1))

  c2 <- Constant(matrix(c(1, 2, 3), ncol = 1))
  expect_true(is_log_log_constant(c2))
})

## @cvxpy NONE
test_that("is_log_log_constant: zero/negative constants are NOT log-log constant", {
  c0 <- Constant(0)
  expect_false(is_log_log_constant(c0))

  cn <- Constant(-1)
  expect_false(is_log_log_constant(cn))

  cm <- Constant(matrix(c(1, -1, 2), ncol = 1))
  expect_false(is_log_log_constant(cm))
})

## @cvxpy NONE
test_that("is_log_log_constant: variables are NOT log-log constant", {
  x <- Variable(3)
  expect_false(is_log_log_constant(x))
})

## @cvxpy NONE
test_that("is_pos on Constant: strictly positive", {
  expect_true(is_pos(Constant(5)))
  expect_true(is_pos(Constant(matrix(c(1, 2, 3), ncol = 1))))
  expect_false(is_pos(Constant(0)))
  expect_false(is_pos(Constant(-1)))
  expect_false(is_pos(Constant(matrix(c(1, 0, 2), ncol = 1))))
})

## @cvxpy NONE
test_that("is_pos on Leaf: uses pos attribute", {
  p <- Variable(3, pos = TRUE)
  expect_true(is_pos(p))

  x <- Variable(3)
  expect_false(is_pos(x))
})

## @cvxpy NONE
test_that("Leaf log-log convex/concave: positive vars", {
  p <- Variable(3, pos = TRUE)
  expect_true(is_log_log_convex(p))
  expect_true(is_log_log_concave(p))
  expect_true(is_log_log_affine(p))

  x <- Variable(3)
  expect_false(is_log_log_convex(x))
  expect_false(is_log_log_concave(x))
  expect_false(is_log_log_affine(x))
})

## @cvxpy NONE
test_that("is_log_log_affine: positive constant is log-log affine", {
  c1 <- Constant(3)
  expect_true(is_log_log_affine(c1))

  c2 <- Constant(0)
  expect_false(is_log_log_affine(c2))
})

## @cvxpy NONE
test_that("log-log convexity of atoms: exp is log-log convex", {
  p <- Variable(3, pos = TRUE)
  e <- exp(p)
  expect_true(is_atom_log_log_convex(e))
  expect_false(is_atom_log_log_concave(e))
})

## @cvxpy NONE
test_that("log-log concavity of atoms: log is log-log concave", {
  p <- Variable(3, pos = TRUE)
  l <- log(p)
  expect_false(is_atom_log_log_convex(l))
  expect_true(is_atom_log_log_concave(l))
})

## @cvxpy NONE
test_that("Affine atoms log-log curvature (CVXPY-aligned)", {
  p <- Variable(3, pos = TRUE)

  ## SumEntries: convex but NOT concave (CVXPY sum.py lines 69-75)
  s <- sum_entries(p)
  expect_true(is_atom_log_log_convex(s))
  expect_false(is_atom_log_log_concave(s))

  ## DiagVec: log-log affine (both convex and concave)
  dv <- DiagVec(p)
  expect_true(is_atom_log_log_convex(dv))
  expect_true(is_atom_log_log_concave(dv))
})

## @cvxpy NONE
test_that("MaxEntries log-log: convex not concave", {
  p <- Variable(3, pos = TRUE)
  m <- max_entries(p)
  expect_true(is_atom_log_log_convex(m))
  expect_false(is_atom_log_log_concave(m))
})

## @cvxpy NONE
test_that("MinEntries log-log: concave not convex", {
  p <- Variable(3, pos = TRUE)
  m <- min_entries(p)
  expect_false(is_atom_log_log_convex(m))
  expect_true(is_atom_log_log_concave(m))
})

## ═══════════════════════════════════════════════════════════════════
## B7: QuadOverLin inherits AxisAtom
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("QuadOverLin: inherits from AxisAtom", {
  x <- Variable(3)
  q <- quad_over_lin(x, Constant(1))
  expect_true(inherits(q, "CVXR::AxisAtom"))
})

## @cvxpy NONE
test_that("QuadOverLin: default axis=NULL gives scalar output", {
  x <- Variable(c(3, 2))
  q <- quad_over_lin(x, Constant(1))
  expect_equal(q@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("QuadOverLin: axis=2 reduces column-wise", {
  x <- Variable(c(3, 2))
  q <- quad_over_lin(x, Constant(1), axis = 2L)
  ## axis=2: sum along rows (column-wise) → (1, ncol) = (1, 2)
  expect_equal(q@shape, c(1L, 2L))
})

## @cvxpy NONE
test_that("QuadOverLin: axis=1 reduces cols", {
  x <- Variable(c(3, 2))
  q <- quad_over_lin(x, Constant(1), axis = 1L)
  ## axis=1: sum along cols -> result has nrow entries -> c(3, 1)
  expect_equal(q@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("QuadOverLin: numeric_value with axis", {
  x <- Variable(c(2, 3))
  xval <- matrix(1:6, nrow = 2, ncol = 3)
  q_null <- quad_over_lin(x, Constant(2))
  q_ax0 <- quad_over_lin(x, Constant(2), axis = 2L)
  q_ax1 <- quad_over_lin(x, Constant(2), axis = 1L)

  ## axis=NULL: sum all squares / y
  nv_null <- numeric_value(q_null, list(xval, matrix(2, 1, 1)))
  expect_equal(nv_null, matrix(sum(xval^2) / 2, 1, 1))

  ## axis=0: column sums of squares / y
  nv_ax0 <- numeric_value(q_ax0, list(xval, matrix(2, 1, 1)))
  expected_ax0 <- apply(xval^2, 2, sum) / 2
  expect_equal(as.numeric(nv_ax0), expected_ax0)

  ## axis=1: row sums of squares / y
  nv_ax1 <- numeric_value(q_ax1, list(xval, matrix(2, 1, 1)))
  expected_ax1 <- apply(xval^2, 1, sum) / 2
  expect_equal(as.numeric(nv_ax1), expected_ax1)
})

## @cvxpy NONE
test_that("QuadOverLin: get_data returns axis and keepdims", {
  x <- Variable(3)
  q <- quad_over_lin(x, Constant(1), axis = 2L, keepdims = TRUE)
  d <- get_data(q)
  expect_equal(d, list(2L, TRUE))
})

## @cvxpy NONE
test_that("QuadOverLin: has_quadratic_term and is_qpwa", {
  x <- Variable(3)
  q <- quad_over_lin(x, Constant(1))
  expect_true(has_quadratic_term(q))
  expect_true(is_qpwa(q))
})

## ═══════════════════════════════════════════════════════════════════
## B8: DiagVec/DiagMat with k parameter
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DiagVec: k=0 (main diagonal) works as before", {
  x <- Variable(3)
  dv <- DiagVec(x)
  expect_equal(dv@shape, c(3L, 3L))
  expect_equal(dv@k, 0L)
  expect_true(is_symmetric(dv))
})

## @cvxpy NONE
test_that("DiagVec: k > 0 (super-diagonal)", {
  x <- Variable(3)
  dv <- DiagVec(x, k = 1L)
  expect_equal(dv@shape, c(4L, 4L))
  expect_equal(dv@k, 1L)
  expect_false(is_symmetric(dv))
  expect_false(is_hermitian(dv))
})

## @cvxpy NONE
test_that("DiagVec: k < 0 (sub-diagonal)", {
  x <- Variable(3)
  dv <- DiagVec(x, k = -2L)
  expect_equal(dv@shape, c(5L, 5L))
  expect_equal(dv@k, -2L)
})

## @cvxpy NONE
test_that("DiagVec: numeric_value with k", {
  x <- Variable(3)
  vals <- matrix(c(1, 2, 3), ncol = 1)

  ## k=0
  dv0 <- DiagVec(x, k = 0L)
  nv0 <- numeric_value(dv0, list(vals))
  expect_equal(diag(nv0), c(1, 2, 3))

  ## k=1 (super-diagonal)
  dv1 <- DiagVec(x, k = 1L)
  nv1 <- numeric_value(dv1, list(vals))
  expect_equal(nv1[1, 2], 1)
  expect_equal(nv1[2, 3], 2)
  expect_equal(nv1[3, 4], 3)
  expect_equal(nv1[1, 1], 0)

  ## k=-1 (sub-diagonal)
  dvm1 <- DiagVec(x, k = -1L)
  nvm1 <- numeric_value(dvm1, list(vals))
  expect_equal(nvm1[2, 1], 1)
  expect_equal(nvm1[3, 2], 2)
  expect_equal(nvm1[4, 3], 3)
})

## @cvxpy NONE
test_that("DiagVec: is_psd and is_nsd with k", {
  x <- Variable(3, nonneg = TRUE)
  dv0 <- DiagVec(x, k = 0L)
  expect_true(is_psd(dv0))

  dv1 <- DiagVec(x, k = 1L)
  expect_false(is_psd(dv1))
})

## @cvxpy NONE
test_that("DiagVec: get_data returns k", {
  x <- Variable(3)
  dv <- DiagVec(x, k = 2L)
  expect_equal(get_data(dv), list(2L))
})

## @cvxpy NONE
test_that("DiagMat: k=0 (main diagonal) works as before", {
  x <- Variable(c(3, 3))
  dm <- DiagMat(x)
  expect_equal(dm@shape, c(3L, 1L))
  expect_equal(dm@k, 0L)
})

## @cvxpy NONE
test_that("DiagMat: k > 0 extracts super-diagonal", {
  x <- Variable(c(4, 4))
  dm <- DiagMat(x, k = 1L)
  expect_equal(dm@shape, c(3L, 1L))
  expect_equal(dm@k, 1L)
})

## @cvxpy NONE
test_that("DiagMat: numeric_value with k", {
  x <- Variable(c(4, 4))
  m <- matrix(1:16, nrow = 4, ncol = 4)

  ## k=0: main diagonal
  dm0 <- DiagMat(x, k = 0L)
  nv0 <- numeric_value(dm0, list(m))
  expect_equal(as.numeric(nv0), c(m[1,1], m[2,2], m[3,3], m[4,4]))

  ## k=1: super-diagonal
  dm1 <- DiagMat(x, k = 1L)
  nv1 <- numeric_value(dm1, list(m))
  expect_equal(as.numeric(nv1), c(m[1,2], m[2,3], m[3,4]))

  ## k=-1: sub-diagonal
  dmn1 <- DiagMat(x, k = -1L)
  nvn1 <- numeric_value(dmn1, list(m))
  expect_equal(as.numeric(nvn1), c(m[2,1], m[3,2], m[4,3]))
})

## @cvxpy NONE
test_that("DiagMat: validation rejects out-of-bounds k", {
  x <- Variable(c(3, 3))
  expect_error(DiagMat(x, k = 3L), "out of bounds")
  expect_error(DiagMat(x, k = -3L), "out of bounds")
})

## @cvxpy NONE
test_that("DiagMat: get_data returns k", {
  x <- Variable(c(4, 4))
  dm <- DiagMat(x, k = -1L)
  expect_equal(get_data(dm), list(-1L))
})

## ═══════════════════════════════════════════════════════════════════
## H1: Trace sign propagation (PSD/NSD)
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Trace: sign from PSD/NSD args", {
  x <- Variable(c(3, 3), PSD = TRUE)
  tr <- matrix_trace(x)
  expect_true(is_nonneg(tr))
})

## ═══════════════════════════════════════════════════════════════════
## H6: UpperTri row-major ordering
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("UpperTri: numeric_value uses row-major ordering", {
  x <- Variable(c(3, 3))
  ## Create matrix where we know the layout:
  ##   1 2 3
  ##   4 5 6
  ##   7 8 9
  m <- matrix(1:9, nrow = 3, byrow = TRUE)
  ut <- upper_tri(x)
  nv <- numeric_value(ut, list(m))
  ## Row-major upper triangle: (1,2)=2, (1,3)=3, (2,3)=6
  expect_equal(as.numeric(nv), c(2, 3, 6))
})

## ═══════════════════════════════════════════════════════════════════
## H14: Huber validation
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Huber: validates M >= 0", {
  x <- Variable(3)
  expect_error(huber(x, M = -1), "non-negative")
})

## @cvxpy NONE
test_that("Huber: default M=1 works", {
  x <- Variable(3)
  h <- huber(x)
  expect_equal(h@shape, c(3L, 1L))
})

## ═══════════════════════════════════════════════════════════════════
## Log-log on QuadOverLin
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("QuadOverLin: log-log convex", {
  x <- Variable(3)
  q <- quad_over_lin(x, Constant(1))
  expect_true(is_atom_log_log_convex(q))
  expect_false(is_atom_log_log_concave(q))
})
