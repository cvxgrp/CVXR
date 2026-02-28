## Phase 3D: Non-affine Atoms Tests

library(testthat)
library(CVXR)

## ── QuadForm ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("QuadForm: basic construction", {
  x <- Variable(3)
  P <- diag(3)
  q <- QuadForm(x, P)
  expect_equal(q@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("QuadForm: numeric value x'Px", {
  x <- Variable(3)
  P <- diag(3)
  value(x) <- c(1, 2, 3)
  q <- QuadForm(x, P)
  v <- value(q)
  expect_equal(as.numeric(v), 14, tolerance = 1e-12)  # 1+4+9=14
})

## @cvxpy NONE
test_that("QuadForm: is_quadratic for affine x", {
  x <- Variable(3)
  P <- diag(3)
  q <- QuadForm(x, P)
  expect_true(is_quadratic(q))
})

## @cvxpy NONE
test_that("quad_form: convenience function", {
  x <- Variable(3)
  P <- diag(3)
  q <- quad_form(x, P)
  expect_true(S7_inherits(q, QuadForm))
})

## ── QuadOverLin ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("QuadOverLin: basic construction", {
  x <- Variable(3)
  y <- Variable(1, nonneg = TRUE)
  q <- QuadOverLin(x, y)
  expect_equal(q@shape, c(1L, 1L))
  expect_true(is_atom_convex(q))
})

## @cvxpy NONE
test_that("QuadOverLin: numeric value", {
  x <- Variable(3)
  y <- Variable(1, nonneg = TRUE)
  value(x) <- c(3, 4, 0)
  value(y) <- 5
  q <- QuadOverLin(x, y)
  v <- value(q)
  expect_equal(as.numeric(v), 5, tolerance = 1e-12)  # (9+16)/5=5
})

## @cvxpy NONE
test_that("QuadOverLin: domain — y >= 0", {
  x <- Variable(3)
  y <- Variable(1)
  q <- QuadOverLin(x, y)
  d <- domain(q)
  expect_true(length(d) >= 1L)
})

## @cvxpy NONE
test_that("quad_over_lin: convenience function", {
  x <- Variable(3)
  y <- Variable(1, nonneg = TRUE)
  q <- quad_over_lin(x, y)
  expect_true(S7_inherits(q, QuadOverLin))
})

## @cvxpy NONE
test_that("sum_squares: convenience = quad_over_lin(x, 1)", {
  x <- Variable(3)
  s <- sum_squares(x)
  value(x) <- c(1, 2, 3)
  expect_equal(as.numeric(value(s)), 14, tolerance = 1e-12)
})

## ── SumLargest ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("SumLargest: basic construction", {
  x <- Variable(5)
  s <- SumLargest(x, k = 3L)
  expect_equal(s@shape, c(1L, 1L))
  expect_true(is_atom_convex(s))
})

## @cvxpy NONE
test_that("SumLargest: numeric value", {
  x <- Variable(5)
  value(x) <- c(1, 5, 2, 4, 3)
  s <- SumLargest(x, k = 2L)
  v <- value(s)
  expect_equal(as.numeric(v), 9)  # 5 + 4 = 9
})

## @cvxpy NONE
test_that("sum_largest: convenience function", {
  x <- Variable(5)
  s <- sum_largest(x, k = 3L)
  expect_true(S7_inherits(s, SumLargest))
})
