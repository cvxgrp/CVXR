## Phase 3C + 3E: Axis Atoms + Affine Atoms Tests

library(testthat)
library(CVXR)

## ── SumEntries ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("SumEntries: axis=NULL reduces to scalar", {
  x <- Variable(c(3, 4))
  s <- SumEntries(x)
  expect_equal(s@shape, c(1L, 1L))
  expect_true(S7_inherits(s, AxisAffAtom))
})

## @cvxpy NONE
test_that("SumEntries: axis=2 sums column-wise", {
  x <- Variable(c(3, 4))
  s <- SumEntries(x, axis = 2L)
  expect_equal(s@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("SumEntries: axis=1 sums over columns", {
  x <- Variable(c(3, 4))
  s <- SumEntries(x, axis = 1L)
  expect_equal(s@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("SumEntries: keepdims=TRUE", {
  x <- Variable(c(3, 4))
  s <- SumEntries(x, axis = 2L, keepdims = TRUE)
  expect_equal(s@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("SumEntries: numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(1:6, 2, 3)
  s <- SumEntries(x)
  expect_equal(value(s), matrix(21, 1, 1))
})

## @cvxpy NONE
test_that("SumEntries: axis=2 numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(1:6, 2, 3)
  s <- SumEntries(x, axis = 2L)
  v <- value(s)
  expect_equal(v, matrix(c(3, 7, 11), nrow = 1))
})

## @cvxpy NONE
test_that("SumEntries: axis=1 numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(1:6, 2, 3)
  s <- SumEntries(x, axis = 1L)
  v <- value(s)
  expect_equal(v, matrix(c(9, 12), ncol = 1))
})

## @cvxpy NONE
test_that("SumEntries: is affine", {
  x <- Variable(c(2, 3))
  s <- SumEntries(x)
  expect_true(is_atom_convex(s))
  expect_true(is_atom_concave(s))
  expect_true(is_affine(s))
})

## @cvxpy NONE
test_that("sum_entries: convenience function", {
  x <- Variable(c(2, 3))
  s <- sum_entries(x, axis = 1L)
  expect_true(S7_inherits(s, SumEntries))
})

## ── MaxEntries ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("MaxEntries: basic construction", {
  x <- Variable(c(3, 4))
  m <- MaxEntries(x)
  expect_equal(m@shape, c(1L, 1L))
  expect_true(is_atom_convex(m))
  expect_false(is_atom_concave(m))
})

## @cvxpy NONE
test_that("MaxEntries: axis=2", {
  x <- Variable(c(3, 4))
  m <- MaxEntries(x, axis = 2L)
  expect_equal(m@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("MaxEntries: numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(c(1, 5, 3, 2, 4, 6), 2, 3)
  m <- MaxEntries(x)
  expect_equal(as.numeric(value(m)), 6)
})

## @cvxpy NONE
test_that("MaxEntries: DCP — MaxEntries(affine) is convex", {
  x <- Variable(c(2, 3))
  m <- MaxEntries(x)
  expect_true(is_dcp(m))
})

## @cvxpy NONE
test_that("max_entries: convenience function", {
  x <- Variable(c(2, 3))
  m <- max_entries(x, axis = 1L)
  expect_true(S7_inherits(m, MaxEntries))
})

## ── MinEntries ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("MinEntries: basic construction", {
  x <- Variable(c(3, 4))
  m <- MinEntries(x)
  expect_equal(m@shape, c(1L, 1L))
  expect_false(is_atom_convex(m))
  expect_true(is_atom_concave(m))
})

## @cvxpy NONE
test_that("MinEntries: numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(c(1, 5, 3, 2, 4, 6), 2, 3)
  m <- MinEntries(x)
  expect_equal(as.numeric(value(m)), 1)
})

## @cvxpy NONE
test_that("MinEntries: DCP — -MinEntries(affine) is convex", {
  x <- Variable(c(2, 3))
  m <- -MinEntries(x)
  expect_true(is_dcp(m))
})

## ── Norm1 ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Norm1: basic construction", {
  x <- Variable(3)
  n <- Norm1(x)
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
  expect_false(is_atom_concave(n))
})

## @cvxpy NONE
test_that("Norm1: with axis=2", {
  x <- Variable(c(3, 4))
  n <- Norm1(x, axis = 2L)
  expect_equal(n@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("Norm1: numeric value", {
  x <- Variable(3)
  value(x) <- c(-1, 2, -3)
  n <- Norm1(x)
  expect_equal(as.numeric(value(n)), 6)
})

## @cvxpy NONE
test_that("norm1: convenience function", {
  x <- Variable(3)
  n <- norm1(x)
  expect_true(S7_inherits(n, Norm1))
})

## ── NormInf ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("NormInf: basic construction", {
  x <- Variable(3)
  n <- NormInf(x)
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
})

## @cvxpy NONE
test_that("NormInf: numeric value", {
  x <- Variable(3)
  value(x) <- c(-1, 2, -3)
  n <- NormInf(x)
  expect_equal(as.numeric(value(n)), 3)
})

## ── Pnorm ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Pnorm: p=2 (Euclidean norm)", {
  x <- Variable(3)
  n <- Pnorm(x, p = 2)
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
})

## @cvxpy NONE
test_that("Pnorm: numeric value p=2", {
  x <- Variable(3)
  value(x) <- c(3, 4, 0)
  n <- Pnorm(x, p = 2)
  expect_equal(as.numeric(value(n)), 5, tolerance = 1e-12)
})

## @cvxpy NONE
test_that("p_norm: convenience function", {
  x <- Variable(3)
  n <- p_norm(x, p = 2)
  expect_true(S7_inherits(n, Pnorm))
})

## ── cvxr_norm: factory function ─────────────────────────────────

## @cvxpy NONE
test_that("cvxr_norm: p=1 gives Norm1", {
  x <- Variable(3)
  n <- cvxr_norm(x, p = 1)
  expect_true(S7_inherits(n, Norm1))
})

## @cvxpy NONE
test_that("cvxr_norm: p=Inf gives NormInf", {
  x <- Variable(3)
  n <- cvxr_norm(x, p = Inf)
  expect_true(S7_inherits(n, NormInf))
})

## @cvxpy NONE
test_that("cvxr_norm: p=2 gives Pnorm", {
  x <- Variable(3)
  n <- cvxr_norm(x, p = 2)
  expect_true(S7_inherits(n, Pnorm))
})

## ── Reshape ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Reshape: basic construction", {
  x <- Variable(c(2, 3))
  r <- Reshape(x, c(6L, 1L))
  expect_equal(r@shape, c(6L, 1L))
  expect_true(S7_inherits(r, AffAtom))
})

## @cvxpy NONE
test_that("Reshape: numeric value", {
  x <- Variable(c(2, 3))
  value(x) <- matrix(1:6, 2, 3)
  r <- Reshape(x, c(3L, 2L))
  v <- value(r)
  expect_equal(dim(v), c(3L, 2L))
})

## @cvxpy NONE
test_that("reshape_expr: convenience function", {
  x <- Variable(c(2, 3))
  r <- reshape_expr(x, c(6L, 1L))
  expect_true(S7_inherits(r, Reshape))
})

## ── DiagVec / DiagMat ───────────────────────────────────────────

## @cvxpy NONE
test_that("DiagVec: vector to diagonal matrix", {
  x <- Variable(3)
  d <- DiagVec(x)
  expect_equal(d@shape, c(3L, 3L))
  expect_true(S7_inherits(d, AffAtom))
})

## @cvxpy NONE
test_that("DiagMat: matrix to diagonal vector", {
  x <- Variable(c(3, 3))
  d <- DiagMat(x)
  expect_equal(d@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("DiagVec: numeric value", {
  x <- Variable(3)
  value(x) <- c(1, 2, 3)
  d <- DiagVec(x)
  v <- value(d)
  expect_equal(diag(v), c(1, 2, 3))
})

## @cvxpy NONE
test_that("DiagMat: numeric value", {
  x <- Variable(c(3, 3))
  value(x) <- diag(c(5, 10, 15))
  d <- DiagMat(x)
  expect_equal(value(d), matrix(c(5, 10, 15), ncol = 1))
})

## ── Trace ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Trace: basic construction", {
  x <- Variable(c(3, 3))
  tr <- Trace(x)
  expect_equal(tr@shape, c(1L, 1L))
  expect_true(S7_inherits(tr, AffAtom))
})

## @cvxpy NONE
test_that("Trace: numeric value", {
  x <- Variable(c(3, 3))
  value(x) <- diag(c(1, 2, 3))
  tr <- Trace(x)
  expect_equal(as.numeric(value(tr)), 6)
})

## @cvxpy NONE
test_that("matrix_trace: convenience function", {
  x <- Variable(c(3, 3))
  tr <- matrix_trace(x)
  expect_true(S7_inherits(tr, Trace))
})

## ── HStack ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HStack: basic construction", {
  x <- Variable(c(3, 2))
  y <- Variable(c(3, 4))
  h <- HStack(x, y)
  expect_equal(h@shape, c(3L, 6L))
  expect_true(S7_inherits(h, AffAtom))
})

## @cvxpy NONE
test_that("HStack: numeric value", {
  x <- Variable(c(2, 2))
  y <- Variable(c(2, 3))
  value(x) <- matrix(1:4, 2, 2)
  value(y) <- matrix(5:10, 2, 3)
  h <- HStack(x, y)
  v <- value(h)
  expect_equal(dim(v), c(2L, 5L))
  expect_equal(v[, 1:2], matrix(1:4, 2, 2))
})

## @cvxpy NONE
test_that("hstack: convenience function", {
  x <- Variable(c(3, 2))
  y <- Variable(c(3, 4))
  h <- hstack(x, y)
  expect_true(S7_inherits(h, HStack))
})

## ── VStack ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("VStack: basic construction", {
  x <- Variable(c(2, 3))
  y <- Variable(c(4, 3))
  v <- VStack(x, y)
  expect_equal(v@shape, c(6L, 3L))
})

## @cvxpy NONE
test_that("VStack: numeric value", {
  x <- Variable(c(2, 2))
  y <- Variable(c(3, 2))
  value(x) <- matrix(1:4, 2, 2)
  value(y) <- matrix(5:10, 3, 2)
  v <- VStack(x, y)
  val <- value(v)
  expect_equal(dim(val), c(5L, 2L))
})

## @cvxpy NONE
test_that("vstack: convenience function", {
  x <- Variable(c(2, 3))
  y <- Variable(c(4, 3))
  v <- vstack(x, y)
  expect_true(S7_inherits(v, VStack))
})

## ── Kron ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Kron: basic construction", {
  A <- matrix(1:4, 2, 2)
  x <- Variable(c(3, 3))
  k <- Kron(A, x)
  expect_equal(k@shape, c(6L, 6L))
})

## @cvxpy NONE
test_that("kron: convenience function", {
  A <- matrix(1:4, 2, 2)
  x <- Variable(c(3, 3))
  k <- kron(A, x)
  expect_true(S7_inherits(k, Kron))
})

## ── UpperTri ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("UpperTri: basic construction", {
  x <- Variable(c(3, 3))
  u <- UpperTri(x)
  ## 3 upper-tri entries (strict)
  expect_equal(u@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("UpperTri: numeric value", {
  x <- Variable(c(3, 3))
  value(x) <- matrix(1:9, 3, 3)
  u <- UpperTri(x)
  v <- value(u)
  ## Strict upper triangle: (1,2)=4, (1,3)=7, (2,3)=8
  expect_equal(length(v), 3L)
})

## @cvxpy NONE
test_that("upper_tri: convenience function", {
  x <- Variable(c(3, 3))
  u <- upper_tri(x)
  expect_true(S7_inherits(u, UpperTri))
})

## ── Convolve ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Convolve: basic construction", {
  a <- c(1, 1, 1)
  x <- Variable(5)
  c_expr <- Convolve(a, x)
  ## Convolution length: 3 + 5 - 1 = 7
  expect_equal(c_expr@shape, c(7L, 1L))
})

## @cvxpy NONE
test_that("conv: convenience function", {
  a <- c(1, 1, 1)
  x <- Variable(5)
  c_expr <- conv(a, x)
  expect_true(S7_inherits(c_expr, Convolve))
})

## ── Cumsum ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Cumsum: basic construction", {
  x <- Variable(5)
  c_atom <- Cumsum(x)
  ## Same shape as input
  expect_equal(c_atom@shape, c(5L, 1L))
  expect_true(S7_inherits(c_atom, AxisAffAtom))
})

## @cvxpy NONE
test_that("Cumsum: numeric value", {
  x <- Variable(4)
  value(x) <- c(1, 2, 3, 4)
  c_atom <- Cumsum(x)
  expect_equal(value(c_atom), matrix(c(1, 3, 6, 10), ncol = 1))
})

## @cvxpy NONE
test_that("cumsum_axis: convenience function", {
  x <- Variable(5)
  c_atom <- cumsum_axis(x)
  expect_true(S7_inherits(c_atom, Cumsum))
})
