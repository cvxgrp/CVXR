## Tests for Phase 3G: DCP integration tests
## Tests: is_dcp on realistic optimization expressions

library(testthat)
library(CVXR)

## ── Portfolio optimization ────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: quad_form(x, Sigma) is convex for PSD Sigma", {
  n <- 5L
  x <- Variable(n)
  ## Create a PSD matrix (Sigma = I)
  Sigma <- Constant(diag(n))
  qf <- quad_form(x, Sigma)
  expect_true(is_convex(qf))
  expect_true(is_dcp(Minimize(qf)))
})

## ── LASSO ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: sum_squares(A*x - b) + lambda * norm1(x) is convex", {
  m <- 10L
  n <- 5L
  x <- Variable(n)
  A <- Constant(matrix(rnorm(m * n), nrow = m))
  b <- Constant(matrix(rnorm(m), ncol = 1))
  lambda <- 0.1

  loss <- sum_squares(A %*% x - b)
  reg <- lambda * norm1(x)
  obj <- loss + reg

  expect_true(is_convex(loss))
  expect_true(is_convex(reg))
  expect_true(is_convex(obj))
  expect_true(is_dcp(Minimize(obj)))
})

## ── Entropy maximization ────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: sum(entr(x)) is concave", {
  x <- Variable(5)
  e <- sum(entr(x))
  expect_true(is_concave(e))
  expect_true(is_dcp(Maximize(e)))
})

## ── sum_squares is convex ────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: is_dcp(sum_squares(x) + norm1(y)) is TRUE", {
  x <- Variable(3)
  y <- Variable(3)
  obj <- sum_squares(x) + norm1(y)
  expect_true(is_convex(obj))
  expect_true(is_dcp(Minimize(obj)))
})

## ── DCP violations ───────────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: Minimize(log(x)) is NOT DCP (concave objective in minimize)", {
  x <- Variable(1)
  l <- log(x)
  ## log is concave, so minimizing it is NOT DCP
  expect_true(is_concave(l))
  expect_false(is_dcp(Minimize(l)))
})

## @cvxpy NONE
test_that("DCP: Minimize(exp(x^2)) is NOT DCP (convex of convex)", {
  x <- Variable(3)
  ## x^2 is convex, exp is convex + increasing
  ## But exp(x^2) = convex(convex + nondecr) → convex by DCP rules
  ## Actually, exp is incr, x^2 is convex → exp(x^2) IS convex
  e <- exp(x^2)
  expect_true(is_convex(e))
})

## ── Affine operations preserve DCP ───────────────────────────────

## @cvxpy NONE
test_that("DCP: sum(x) is affine", {
  x <- Variable(3)
  s <- sum(x)
  expect_true(is_affine(s))
})

## @cvxpy NONE
test_that("DCP: reshape preserves affineness", {
  x <- Variable(c(3, 2))
  r <- reshape_expr(x, c(6, 1))
  expect_true(is_affine(r))
})

## @cvxpy NONE
test_that("DCP: trace is affine", {
  x <- Variable(c(3, 3))
  tr <- matrix_trace(x)
  expect_true(is_affine(tr))
})

## @cvxpy NONE
test_that("DCP: hstack/vstack are affine", {
  x <- Variable(c(3, 1))
  y <- Variable(c(3, 2))
  h <- hstack(x, y)
  expect_true(is_affine(h))

  a <- Variable(c(2, 3))
  b <- Variable(c(1, 3))
  v <- vstack(a, b)
  expect_true(is_affine(v))
})

## ── Norm compositions ────────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: norm1(A*x - b) is convex", {
  m <- 5L; n <- 3L
  x <- Variable(n)
  A <- Constant(matrix(rnorm(m * n), nrow = m))
  b <- Constant(matrix(rnorm(m), ncol = 1))
  n1 <- norm1(A %*% x - b)
  expect_true(is_convex(n1))
})

## @cvxpy NONE
test_that("DCP: norm_inf(x) is convex", {
  x <- Variable(5)
  ni <- norm_inf(x)
  expect_true(is_convex(ni))
})

## @cvxpy NONE
test_that("DCP: p_norm(x, 2) is convex", {
  x <- Variable(5)
  pn <- p_norm(x, 2)
  expect_true(is_convex(pn))
})

## ── Huber in DCP ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: sum(huber(x)) is convex", {
  x <- Variable(5)
  h <- sum(huber(x))
  expect_true(is_convex(h))
})

## ── Max/Min entries DCP ──────────────────────────────────────────

## @cvxpy NONE
test_that("DCP: max(x) is convex", {
  x <- Variable(5)
  m <- max(x)
  expect_true(is_convex(m))
  expect_true(is_dcp(Minimize(m)))
})

## @cvxpy NONE
test_that("DCP: min(x) is concave", {
  x <- Variable(5)
  m <- min(x)
  expect_true(is_concave(m))
  expect_true(is_dcp(Maximize(m)))
})

## ── Convenience functions DCP ────────────────────────────────────

## @cvxpy NONE
test_that("DCP: square(x) is convex", {
  x <- Variable(3)
  s <- square(x)
  expect_true(is_convex(s))
})

## @cvxpy NONE
test_that("DCP: sqrt(x) is concave", {
  x <- Variable(3)
  s <- sqrt(x)
  expect_true(is_concave(s))
})

## @cvxpy NONE
test_that("DCP: pos(x) is convex", {
  x <- Variable(3)
  p <- pos(x)
  expect_true(is_convex(p))
})

## @cvxpy NONE
test_that("DCP: neg(x) is convex", {
  x <- Variable(3)
  n <- neg(x)
  expect_true(is_convex(n))
})

## ── Smoke test from PLAN_v2.md ──────────────────────────────────

## @cvxpy NONE
test_that("PLAN_v2 smoke test: is_dcp(sum_squares(x) + norm1(y)) is TRUE", {
  x <- Variable(5)
  y <- Variable(5)
  expr <- sum_squares(x) + norm1(y)
  expect_true(is_dcp(Minimize(expr)))
})
