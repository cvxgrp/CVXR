## Tests for Phase 3C: AxisAtom-based atoms
## Tests: Norm1, NormInf, Pnorm, MaxEntries, MinEntries, cvxr_norm

library(testthat)

## ── Norm1 ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Norm1: scalar output (axis=NULL)", {
  x <- Variable(c(3, 2))
  n <- Norm1(x)
  expect_s3_class(n, "CVXR::Norm1")
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
  expect_true(is_nonneg(n))
})

## @cvxpy NONE
test_that("Norm1: axis=2 reduces column-wise", {
  x <- Variable(c(3, 2))
  n <- Norm1(x, axis = 2L)
  ## axis=2: reduce rows (column-wise) → (1, ncol) = (1, 2)
  expect_equal(n@shape, c(1L, 2L))
})

## @cvxpy NONE
test_that("Norm1: axis=1 reduces columns", {
  x <- Variable(c(3, 2))
  n <- Norm1(x, axis = 1L)
  ## axis=1: reduce cols → (nrow, 1) = (3, 1)
  expect_equal(n@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Norm1: numeric evaluation", {
  x <- Variable(3)
  n <- Norm1(x)
  val <- matrix(c(-2, 3, -1), ncol = 1)
  result <- numeric_value(n, list(val))
  expect_equal(as.numeric(result), 6)
})

## @cvxpy NONE
test_that("norm1() factory function", {
  x <- Variable(3)
  n <- norm1(x)
  expect_s3_class(n, "CVXR::Norm1")
})

## ── NormInf ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("NormInf: basic properties", {
  x <- Variable(3)
  n <- NormInf(x)
  expect_s3_class(n, "CVXR::NormInf")
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
  expect_true(is_nonneg(n))
})

## @cvxpy NONE
test_that("NormInf: numeric evaluation", {
  x <- Variable(3)
  n <- NormInf(x)
  val <- matrix(c(-5, 3, -1), ncol = 1)
  result <- numeric_value(n, list(val))
  expect_equal(as.numeric(result), 5)
})

## @cvxpy NONE
test_that("norm_inf() factory function", {
  x <- Variable(3)
  n <- norm_inf(x)
  expect_s3_class(n, "CVXR::NormInf")
})

## ── Pnorm ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Pnorm: p=2 (Euclidean)", {
  x <- Variable(3)
  n <- Pnorm(x, 2)
  expect_s3_class(n, "CVXR::Pnorm")
  expect_equal(n@shape, c(1L, 1L))
  expect_true(is_atom_convex(n))
  expect_true(is_nonneg(n))
})

## @cvxpy NONE
test_that("Pnorm: p=2 numeric", {
  x <- Variable(3)
  n <- Pnorm(x, 2)
  val <- matrix(c(3, 4, 0), ncol = 1)
  result <- numeric_value(n, list(val))
  expect_equal(as.numeric(result), 5)
})

## @cvxpy NONE
test_that("p_norm() factory function", {
  x <- Variable(3)
  n <- p_norm(x, 2)
  expect_s3_class(n, "CVXR::Pnorm")
})

## ── MaxEntries ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("MaxEntries: scalar output (axis=NULL)", {
  x <- Variable(c(3, 2))
  m <- MaxEntries(x)
  expect_s3_class(m, "CVXR::MaxEntries")
  expect_equal(m@shape, c(1L, 1L))
  expect_true(is_atom_convex(m))
  expect_false(is_atom_concave(m))
})

## @cvxpy NONE
test_that("MaxEntries: always increasing", {
  x <- Variable(3)
  m <- MaxEntries(x)
  expect_true(is_incr(m, 1L))
  expect_false(is_decr(m, 1L))
})

## @cvxpy NONE
test_that("MaxEntries: numeric evaluation", {
  x <- Variable(3)
  m <- MaxEntries(x)
  val <- matrix(c(-2, 7, 3), ncol = 1)
  result <- numeric_value(m, list(val))
  expect_equal(as.numeric(result), 7)
})

## @cvxpy NONE
test_that("MaxEntries: axis=2", {
  x <- Variable(c(3, 2))
  m <- MaxEntries(x, axis = 2L)
  expect_equal(m@shape, c(1L, 2L))
  val <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  result <- numeric_value(m, list(val))
  expect_equal(as.numeric(result), c(3, 6))
})

## @cvxpy NONE
test_that("max_entries() factory function", {
  x <- Variable(3)
  m <- max_entries(x)
  expect_s3_class(m, "CVXR::MaxEntries")
})

## ── MinEntries ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("MinEntries: scalar output (axis=NULL)", {
  x <- Variable(c(3, 2))
  m <- MinEntries(x)
  expect_s3_class(m, "CVXR::MinEntries")
  expect_equal(m@shape, c(1L, 1L))
  expect_false(is_atom_convex(m))
  expect_true(is_atom_concave(m))
})

## @cvxpy NONE
test_that("MinEntries: numeric evaluation", {
  x <- Variable(3)
  m <- MinEntries(x)
  val <- matrix(c(-2, 7, 3), ncol = 1)
  result <- numeric_value(m, list(val))
  expect_equal(as.numeric(result), -2)
})

## @cvxpy NONE
test_that("min_entries() factory function", {
  x <- Variable(3)
  m <- min_entries(x)
  expect_s3_class(m, "CVXR::MinEntries")
})

## ── cvxr_norm factory ──────────────────────────────────────────────

## @cvxpy NONE
test_that("cvxr_norm dispatches to correct atom", {
  x <- Variable(3)
  ## p=1 → Norm1
  n1 <- cvxr_norm(x, 1)
  expect_s3_class(n1, "CVXR::Norm1")
  ## p=2 → Pnorm
  n2 <- cvxr_norm(x, 2)
  expect_s3_class(n2, "CVXR::Pnorm")
  ## p=Inf → NormInf
  ninf <- cvxr_norm(x, Inf)
  expect_s3_class(ninf, "CVXR::NormInf")
})
