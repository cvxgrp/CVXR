## Tests for Phase 3E: Affine atoms
## Tests: SumEntries, Reshape, DiagVec/DiagMat, Trace, HStack, VStack,
##        Kron, UpperTri, Convolve, Cumsum

library(testthat)

## ── SumEntries ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("SumEntries: axis=NULL → scalar", {
  x <- Variable(c(3, 2))
  s <- SumEntries(x)
  expect_s3_class(s, "CVXR::SumEntries")
  expect_equal(s@shape, c(1L, 1L))
  ## AffAtom: both convex and concave
  expect_true(is_convex(s))
  expect_true(is_concave(s))
})

## @cvxpy NONE
test_that("SumEntries: axis=2 reduces column-wise", {
  x <- Variable(c(3, 2))
  s <- SumEntries(x, axis = 2L)
  expect_equal(s@shape, c(1L, 2L))
})

## @cvxpy NONE
test_that("SumEntries: axis=1 reduces columns", {
  x <- Variable(c(3, 2))
  s <- SumEntries(x, axis = 1L)
  expect_equal(s@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("SumEntries: numeric evaluation", {
  x <- Variable(c(2, 3))
  s <- SumEntries(x)
  val <- matrix(1:6, nrow = 2, ncol = 3)
  result <- numeric_value(s, list(val))
  expect_equal(as.numeric(result), 21)
})

## @cvxpy NONE
test_that("SumEntries: numeric with axis=0", {
  x <- Variable(c(2, 3))
  s <- SumEntries(x, axis = 2L)
  val <- matrix(1:6, nrow = 2, ncol = 3)
  result <- numeric_value(s, list(val))
  expect_equal(as.numeric(result), c(3, 7, 11))
})

## @cvxpy NONE
test_that("sum_entries() factory function", {
  x <- Variable(3)
  s <- sum_entries(x)
  expect_s3_class(s, "CVXR::SumEntries")
})

## ── Reshape ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Reshape: basic construction", {
  x <- Variable(c(3, 2))
  r <- Reshape(x, c(6, 1))
  expect_s3_class(r, "CVXR::Reshape")
  expect_equal(r@shape, c(6L, 1L))
  expect_true(is_convex(r))
  expect_true(is_concave(r))
})

## @cvxpy NONE
test_that("Reshape: -1 dimension inference", {
  x <- Variable(c(3, 4))
  r <- Reshape(x, c(-1, 6))
  expect_equal(r@shape, c(2L, 6L))
})

## @cvxpy NONE
test_that("Reshape: numeric evaluation (column-major)", {
  x <- Variable(c(2, 3))
  r <- Reshape(x, c(3, 2))
  val <- matrix(1:6, nrow = 2, ncol = 3)
  result <- numeric_value(r, list(val))
  expect_equal(dim(result), c(3, 2))
  ## Column-major: as.vector(val) → c(1,2,3,4,5,6) → matrix 3x2
  expect_equal(result[1, 1], 1)
  expect_equal(result[3, 2], 6)
})

## @cvxpy NONE
test_that("Reshape: size mismatch errors", {
  x <- Variable(c(3, 2))
  expect_error(Reshape(x, c(4, 2)), "does not match")
})

## @cvxpy NONE
test_that("reshape_expr() factory function", {
  x <- Variable(c(3, 2))
  r <- reshape_expr(x, c(2, 3))
  expect_s3_class(r, "CVXR::Reshape")
})

## ── DiagVec ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("DiagVec: vector to diagonal matrix", {
  x <- Variable(3)
  d <- DiagVec(x)
  expect_s3_class(d, "CVXR::DiagVec")
  expect_equal(d@shape, c(3L, 3L))
  expect_true(is_symmetric(d))
})

## @cvxpy NONE
test_that("DiagVec: numeric evaluation", {
  x <- Variable(3)
  d <- DiagVec(x)
  val <- matrix(c(1, 2, 3), ncol = 1)
  result <- numeric_value(d, list(val))
  expect_equal(result, diag(c(1, 2, 3)))
})

## @cvxpy NONE
test_that("DiagVec: requires column vector", {
  x <- Variable(c(3, 2))
  expect_error(DiagVec(x), "column vector")
})

## ── DiagMat ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("DiagMat: matrix diagonal to vector", {
  x <- Variable(c(3, 3))
  d <- DiagMat(x)
  expect_s3_class(d, "CVXR::DiagMat")
  expect_equal(d@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("DiagMat: numeric evaluation", {
  x <- Variable(c(3, 3))
  d <- DiagMat(x)
  val <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), 3, 3)
  result <- numeric_value(d, list(val))
  expect_equal(as.numeric(result), c(1, 2, 3))
})

## @cvxpy NONE
test_that("DiagMat: requires square matrix", {
  x <- Variable(c(3, 2))
  expect_error(DiagMat(x), "square matrix")
})

## ── Trace ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Trace: basic construction", {
  x <- Variable(c(3, 3))
  tr <- Trace(x)
  expect_s3_class(tr, "CVXR::Trace")
  expect_equal(tr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Trace: numeric evaluation", {
  x <- Variable(c(3, 3))
  tr <- Trace(x)
  val <- diag(c(2, 3, 5))
  result <- numeric_value(tr, list(val))
  expect_equal(as.numeric(result), 10)
})

## @cvxpy NONE
test_that("matrix_trace() factory function", {
  x <- Variable(c(3, 3))
  tr <- matrix_trace(x)
  expect_s3_class(tr, "CVXR::Trace")
})

## ── HStack ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HStack: horizontal concatenation", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 2))
  h <- HStack(x, y)
  expect_s3_class(h, "CVXR::HStack")
  expect_equal(h@shape, c(3L, 3L))
})

## @cvxpy NONE
test_that("HStack: numeric evaluation", {
  x <- Variable(c(2, 1))
  y <- Variable(c(2, 1))
  h <- HStack(x, y)
  v1 <- matrix(c(1, 2), ncol = 1)
  v2 <- matrix(c(3, 4), ncol = 1)
  result <- numeric_value(h, list(v1, v2))
  expect_equal(result, matrix(c(1, 2, 3, 4), nrow = 2))
})

## @cvxpy NONE
test_that("HStack: row mismatch errors", {
  x <- Variable(c(2, 1))
  y <- Variable(c(3, 1))
  expect_error(HStack(x, y), "row")
})

## @cvxpy NONE
test_that("hstack() factory function", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 2))
  h <- hstack(x, y)
  expect_s3_class(h, "CVXR::HStack")
})

## ── VStack ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("VStack: vertical concatenation", {
  x <- Variable(c(2, 3))
  y <- Variable(c(1, 3))
  v <- VStack(x, y)
  expect_s3_class(v, "CVXR::VStack")
  expect_equal(v@shape, c(3L, 3L))
})

## @cvxpy NONE
test_that("VStack: numeric evaluation", {
  x <- Variable(c(1, 2))
  y <- Variable(c(1, 2))
  v <- VStack(x, y)
  v1 <- matrix(c(1, 2), nrow = 1)
  v2 <- matrix(c(3, 4), nrow = 1)
  result <- numeric_value(v, list(v1, v2))
  expect_equal(result, matrix(c(1, 3, 2, 4), nrow = 2))
})

## @cvxpy NONE
test_that("vstack() factory function", {
  x <- Variable(c(2, 3))
  y <- Variable(c(1, 3))
  v <- vstack(x, y)
  expect_s3_class(v, "CVXR::VStack")
})

## ── Kron ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Kron: basic construction", {
  a <- Constant(diag(2))
  x <- Variable(c(3, 3))
  k <- Kron(a, x)
  expect_s3_class(k, "CVXR::Kron")
  expect_equal(k@shape, c(6L, 6L))
})

## @cvxpy NONE
test_that("Kron: numeric evaluation", {
  a <- Constant(diag(2))
  x <- Variable(c(2, 2))
  k <- Kron(a, x)
  av <- diag(2)
  xv <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- numeric_value(k, list(av, xv))
  expected <- kronecker(av, xv)
  expect_equal(result, expected)
})

## @cvxpy NONE
test_that("Kron: requires one constant", {
  x <- Variable(c(2, 2))
  y <- Variable(c(2, 2))
  expect_error(Kron(x, y), "constant")
})

## @cvxpy NONE
test_that("kron() factory function", {
  a <- Constant(diag(2))
  x <- Variable(c(3, 3))
  k <- kron(a, x)
  expect_s3_class(k, "CVXR::Kron")
})

## ── UpperTri ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("UpperTri: basic construction", {
  x <- Variable(c(3, 3))
  u <- UpperTri(x)
  expect_s3_class(u, "CVXR::UpperTri")
  ## Strict upper triangle of 3x3: 3 elements
  expect_equal(u@shape[1L], 3L)
  expect_equal(u@shape[2L], 1L)
})

## @cvxpy NONE
test_that("UpperTri: numeric evaluation", {
  x <- Variable(c(3, 3))
  u <- UpperTri(x)
  val <- matrix(1:9, 3, 3)
  result <- numeric_value(u, list(val))
  ## Upper triangle (column-major): (1,2), (1,3), (2,3) → vals at (4, 7, 8)
  expect_equal(as.numeric(result), c(4, 7, 8))
})

## @cvxpy NONE
test_that("upper_tri() factory function", {
  x <- Variable(c(3, 3))
  u <- upper_tri(x)
  expect_s3_class(u, "CVXR::UpperTri")
})

## ── Convolve ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Convolve: basic construction", {
  a <- Constant(matrix(c(1, 2, 3), ncol = 1))
  x <- Variable(4)
  cv <- Convolve(a, x)
  expect_s3_class(cv, "CVXR::Convolve")
  ## Output length = 3 + 4 - 1 = 6
  expect_equal(cv@shape, c(6L, 1L))
})

## @cvxpy NONE
test_that("Convolve: requires one constant", {
  x <- Variable(3)
  y <- Variable(4)
  expect_error(Convolve(x, y), "constant")
})

## @cvxpy NONE
test_that("Convolve: numeric evaluation", {
  a <- Constant(matrix(c(1, 1), ncol = 1))
  x <- Variable(3)
  cv <- Convolve(a, x)
  av <- matrix(c(1, 1), ncol = 1)
  xv <- matrix(c(1, 2, 3), ncol = 1)
  result <- numeric_value(cv, list(av, xv))
  ## convolve(c(1,1), c(1,2,3), type="open") = c(1, 3, 5, 3)
  expect_equal(as.numeric(result), c(1, 3, 5, 3))
})

## @cvxpy NONE
test_that("conv() factory function", {
  a <- Constant(matrix(c(1, 2), ncol = 1))
  x <- Variable(3)
  cv <- conv(a, x)
  expect_s3_class(cv, "CVXR::Convolve")
})

## ── Cumsum ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Cumsum: preserves shape", {
  x <- Variable(c(3, 2))
  cs <- Cumsum(x)
  expect_s3_class(cs, "CVXR::Cumsum")
  expect_equal(cs@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Cumsum: numeric evaluation (axis=NULL)", {
  x <- Variable(c(2, 2))
  cs <- Cumsum(x)
  val <- matrix(c(1, 2, 3, 4), nrow = 2)
  result <- numeric_value(cs, list(val))
  ## Column-major: cumsum(c(1,2,3,4)) = c(1,3,6,10) → 2x2
  expect_equal(as.numeric(result), c(1, 3, 6, 10))
})

## @cvxpy NONE
test_that("cumsum_axis() factory function", {
  x <- Variable(3)
  cs <- cumsum_axis(x)
  expect_s3_class(cs, "CVXR::Cumsum")
})
