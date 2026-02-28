## Tests for Phase 3B: Elementwise atoms + power_tools
## Tests: Abs, Exp, Log, Power, Entr, Huber, Maximum, Minimum, convenience

library(testthat)

## ── Abs ────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Abs: basic construction and properties", {
  x <- Variable(3)
  a <- Abs(x)
  expect_s3_class(a, "CVXR::Abs")
  expect_equal(a@shape, c(3L, 1L))
  expect_true(is_atom_convex(a))
  expect_false(is_atom_concave(a))
})

## @cvxpy NONE
test_that("Abs: sign is always nonneg", {
  x <- Variable(3)
  a <- Abs(x)
  expect_true(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("Abs: monotonicity depends on arg sign", {
  x <- Variable(3)
  a <- Abs(x)
  ## For general variable: neither increasing nor decreasing
  expect_false(is_incr(a, 1L))
  expect_false(is_decr(a, 1L))
})

## @cvxpy NONE
test_that("Abs: numeric evaluation", {
  x <- Variable(3)
  a <- Abs(x)
  val <- matrix(c(-2, 3, -1), ncol = 1)
  result <- numeric_value(a, list(val))
  expect_equal(result, matrix(c(2, 3, 1), ncol = 1))
})

## @cvxpy NONE
test_that("Abs: is_pwl when arg is PWL", {
  x <- Variable(3)
  a <- Abs(x)
  ## Variable is affine, hence PWL
  expect_true(is_pwl(a))
})

## @cvxpy NONE
test_that("Abs: convex of affine => DCP compliant", {
  x <- Variable(3)
  a <- Abs(x)
  expect_true(is_convex(a))
})

## ── Exp ────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Exp: basic construction and properties", {
  x <- Variable(c(2, 3))
  e <- Exp(x)
  expect_s3_class(e, "CVXR::Exp")
  expect_equal(e@shape, c(2L, 3L))
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))
})

## @cvxpy NONE
test_that("Exp: sign is always positive", {
  x <- Variable(3)
  e <- Exp(x)
  expect_true(is_nonneg(e))
  expect_false(is_nonpos(e))
})

## @cvxpy NONE
test_that("Exp: always increasing", {
  x <- Variable(3)
  e <- Exp(x)
  expect_true(is_incr(e, 1L))
  expect_false(is_decr(e, 1L))
})

## @cvxpy NONE
test_that("Exp: numeric evaluation", {
  x <- Variable(2)
  e <- Exp(x)
  val <- matrix(c(0, 1), ncol = 1)
  result <- numeric_value(e, list(val))
  expect_equal(result, matrix(c(1, exp(1)), ncol = 1))
})

## @cvxpy NONE
test_that("Exp: convex + increasing arg => DCP compliant for affine arg", {
  x <- Variable(3)
  e <- Exp(x)
  expect_true(is_convex(e))
})

## ── Log ────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Log: basic construction and properties", {
  x <- Variable(3)
  l <- Log(x)
  expect_s3_class(l, "CVXR::Log")
  expect_equal(l@shape, c(3L, 1L))
  expect_false(is_atom_convex(l))
  expect_true(is_atom_concave(l))
})

## @cvxpy NONE
test_that("Log: sign is unknown", {
  x <- Variable(3)
  l <- Log(x)
  expect_false(is_nonneg(l))
  expect_false(is_nonpos(l))
})

## @cvxpy NONE
test_that("Log: domain includes arg >= 0", {
  x <- Variable(3)
  l <- Log(x)
  d <- domain(l)
  expect_length(d, 1L)
  expect_s3_class(d[[1L]], "CVXR::Constraint")
})

## @cvxpy NONE
test_that("Log: numeric evaluation", {
  x <- Variable(2)
  l <- Log(x)
  val <- matrix(c(1, exp(1)), ncol = 1)
  result <- numeric_value(l, list(val))
  expect_equal(result, matrix(c(0, 1), ncol = 1))
})

## @cvxpy NONE
test_that("Log: concave + increasing => DCP concave for affine arg", {
  x <- Variable(3)
  l <- Log(x)
  expect_true(is_concave(l))
})

## ── Power ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Power: x^2 is convex", {
  x <- Variable(3)
  p <- Power(x, 2)
  expect_s3_class(p, "CVXR::Power")
  expect_equal(p@shape, c(3L, 1L))
  expect_true(is_atom_convex(p))
  expect_false(is_atom_concave(p))
})

## @cvxpy NONE
test_that("Power: x^0.5 is concave", {
  x <- Variable(3)
  p <- Power(x, 0.5)
  expect_false(is_atom_convex(p))
  expect_true(is_atom_concave(p))
})

## @cvxpy NONE
test_that("Power: x^1 is affine", {
  x <- Variable(3)
  p <- Power(x, 1)
  expect_true(is_atom_convex(p))
  expect_true(is_atom_concave(p))
})

## @cvxpy NONE
test_that("Power: x^0 is constant", {
  x <- Variable(3)
  p <- Power(x, 0)
  expect_true(is_constant(p))
})

## @cvxpy NONE
test_that("Power: x^(-1) is convex", {
  x <- Variable(3)
  p <- Power(x, -1)
  expect_true(is_atom_convex(p))
  expect_false(is_atom_concave(p))
})

## @cvxpy NONE
test_that("Power: sign", {
  x <- Variable(3)
  ## x^2 is always nonneg
  p2 <- Power(x, 2)
  expect_true(is_nonneg(p2))
  ## x^1 has same sign as x
  p1 <- Power(x, 1)
  expect_false(is_nonneg(p1))  # unknown for general variable
})

## @cvxpy NONE
test_that("Power: numeric evaluation", {
  x <- Variable(2)
  p <- Power(x, 2)
  val <- matrix(c(3, -2), ncol = 1)
  result <- numeric_value(p, list(val))
  expect_equal(result, matrix(c(9, 4), ncol = 1))
})

## @cvxpy NONE
test_that("Power: domain for non-power-of-2 exponent", {
  x <- Variable(3)
  p <- Power(x, 3)  # 3 is not power of 2 → needs x >= 0
  d <- domain(p)
  expect_length(d, 1L)
  expect_s3_class(d[[1L]], "CVXR::Constraint")
})

## @cvxpy NONE
test_that("Power: domain for power-of-2 exponent", {
  x <- Variable(3)
  p <- Power(x, 2)  # 2 is power of 2 → no domain constraint
  d <- atom_domain(p)
  expect_length(d, 0L)
})

## @cvxpy NONE
test_that("Power: is_quadratic", {
  x <- Variable(3)
  ## x^2 where x is affine → is_quadratic
  p <- Power(x, 2)
  expect_true(is_quadratic(p))
  ## x^3 where x is affine → not quadratic (unless constant)
  p3 <- Power(x, 3)
  expect_false(is_quadratic(p3))
})

## @cvxpy NONE
test_that("power() factory function", {
  x <- Variable(3)
  p <- power(x, 2)
  expect_s3_class(p, "CVXR::Power")
})

## ── Entr ────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Entr: basic construction and properties", {
  x <- Variable(3)
  e <- Entr(x)
  expect_s3_class(e, "CVXR::Entr")
  expect_equal(e@shape, c(3L, 1L))
  expect_false(is_atom_convex(e))
  expect_true(is_atom_concave(e))
})

## @cvxpy NONE
test_that("Entr: domain includes arg >= 0", {
  x <- Variable(3)
  e <- Entr(x)
  d <- domain(e)
  expect_length(d, 1L)
})

## @cvxpy NONE
test_that("Entr: numeric evaluation", {
  x <- Variable(2)
  e <- Entr(x)
  val <- matrix(c(1, exp(1)), ncol = 1)
  result <- numeric_value(e, list(val))
  ## -x*log(x): -1*0 = 0, -e*1 = -e
  expect_equal(result[1, 1], 0, tolerance = 1e-10)
  expect_equal(result[2, 1], -exp(1), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("entr() factory function", {
  x <- Variable(3)
  e <- entr(x)
  expect_s3_class(e, "CVXR::Entr")
})

## ── Huber ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Huber: basic construction and properties", {
  x <- Variable(3)
  h <- Huber(x)
  expect_s3_class(h, "CVXR::Huber")
  expect_equal(h@shape, c(3L, 1L))
  expect_true(is_atom_convex(h))
  expect_false(is_atom_concave(h))
})

## @cvxpy NONE
test_that("Huber: sign is always nonneg", {
  x <- Variable(3)
  h <- Huber(x)
  expect_true(is_nonneg(h))
  expect_false(is_nonpos(h))
})

## @cvxpy NONE
test_that("Huber: numeric evaluation with M=1", {
  x <- Variable(3)
  h <- Huber(x, 1)
  val <- matrix(c(0.5, 2, -1), ncol = 1)
  result <- numeric_value(h, list(val))
  ## |0.5| <= 1: 0.5^2 = 0.25
  ## |2| > 1: 2*1*2 - 1 = 3
  ## |-1| <= 1: 1
  expect_equal(as.numeric(result), c(0.25, 3.0, 1.0))
})

## @cvxpy NONE
test_that("huber() factory function", {
  x <- Variable(3)
  h <- huber(x, 2)
  expect_s3_class(h, "CVXR::Huber")
})

## ── Maximum ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Maximum: basic construction", {
  x <- Variable(3)
  y <- Variable(3)
  m <- Maximum(x, y)
  expect_s3_class(m, "CVXR::Maximum")
  expect_equal(m@shape, c(3L, 1L))
  expect_true(is_atom_convex(m))
  expect_false(is_atom_concave(m))
})

## @cvxpy NONE
test_that("Maximum: requires at least 2 args", {
  x <- Variable(3)
  expect_error(Maximum(x), "at least 2")
})

## @cvxpy NONE
test_that("Maximum: numeric evaluation", {
  x <- Variable(2)
  y <- Variable(2)
  m <- Maximum(x, y)
  v1 <- matrix(c(1, 4), ncol = 1)
  v2 <- matrix(c(3, 2), ncol = 1)
  result <- numeric_value(m, list(v1, v2))
  expect_equal(as.numeric(result), c(3, 4))
})

## @cvxpy NONE
test_that("max_elemwise() factory function", {
  x <- Variable(3)
  y <- Variable(3)
  m <- max_elemwise(x, y)
  expect_s3_class(m, "CVXR::Maximum")
})

## ── Minimum ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Minimum: basic construction", {
  x <- Variable(3)
  y <- Variable(3)
  m <- Minimum(x, y)
  expect_s3_class(m, "CVXR::Minimum")
  expect_equal(m@shape, c(3L, 1L))
  expect_false(is_atom_convex(m))
  expect_true(is_atom_concave(m))
})

## @cvxpy NONE
test_that("Minimum: numeric evaluation", {
  x <- Variable(2)
  y <- Variable(2)
  m <- Minimum(x, y)
  v1 <- matrix(c(1, 4), ncol = 1)
  v2 <- matrix(c(3, 2), ncol = 1)
  result <- numeric_value(m, list(v1, v2))
  expect_equal(as.numeric(result), c(1, 2))
})

## @cvxpy NONE
test_that("min_elemwise() factory function", {
  x <- Variable(3)
  y <- Variable(3)
  m <- min_elemwise(x, y)
  expect_s3_class(m, "CVXR::Minimum")
})

## ── Convenience functions ──────────────────────────────────────────

## @cvxpy NONE
test_that("square(x) creates Power with p=2", {
  x <- Variable(3)
  s <- square(x)
  expect_s3_class(s, "CVXR::Power")
  expect_equal(s@p_used, 2)
})

## @cvxpy NONE
test_that("pos(x) creates Maximum(x, 0)", {
  x <- Variable(3)
  p <- pos(x)
  expect_s3_class(p, "CVXR::Maximum")
})

## @cvxpy NONE
test_that("neg(x) creates negation of Minimum(x, 0)", {
  x <- Variable(3)
  n <- neg(x)
  ## neg(x) = -min(x, 0) = NegExpression(Minimum(x, 0))
  expect_s3_class(n, "CVXR::NegExpression")
})

## @cvxpy NONE
test_that("inv_pos(x) creates Power with p=-1", {
  x <- Variable(3)
  ip <- inv_pos(x)
  expect_s3_class(ip, "CVXR::Power")
  expect_equal(ip@p_used, -1)
})

## ── power_tools ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("is_power2 correctly identifies powers of 2", {
  expect_true(is_power2(1))
  expect_true(is_power2(2))
  expect_true(is_power2(4))
  expect_true(is_power2(8))
  expect_false(is_power2(3))
  expect_false(is_power2(5))
  expect_false(is_power2(0))
  expect_false(is_power2(-1))
  expect_false(is_power2(1.5))
})
