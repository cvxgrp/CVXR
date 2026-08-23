## Tests for Phase 3F: S3 Math/Summary group dispatch
## Tests: abs(), exp(), log(), sqrt(), cumsum(), sum(), max(), min(), x^p

library(testthat)

## ── Math group: abs ──────────────────────────────────────────────

## @cvxpy NONE
test_that("abs(x) creates Abs atom via Math dispatch", {
  x <- Variable(3)
  a <- abs(x)
  expect_s3_class(a, "CVXR::Abs")
  expect_equal(a@shape, c(3L, 1L))
  expect_true(is_convex(a))
})

## ── Math group: exp ──────────────────────────────────────────────

## @cvxpy NONE
test_that("exp(x) creates Exp atom via Math dispatch", {
  x <- Variable(c(2, 3))
  e <- exp(x)
  expect_s3_class(e, "CVXR::Exp")
  expect_equal(e@shape, c(2L, 3L))
  expect_true(is_convex(e))
})

## ── Math group: log ──────────────────────────────────────────────

## @cvxpy NONE
test_that("log(x) creates Log atom via Math dispatch", {
  x <- Variable(3)
  l <- log(x)
  expect_s3_class(l, "CVXR::Log")
  expect_equal(l@shape, c(3L, 1L))
  expect_true(is_concave(l))
})

## ── Math group: sqrt ─────────────────────────────────────────────

## @cvxpy NONE
test_that("sqrt(x) creates PowerApprox(x, 0.5) via Math dispatch", {
  x <- Variable(3)
  s <- sqrt(x)
  ## sqrt.py:19 is `power(x, Fraction(1, 2))` -- the WRAPPER, whose
  ## `approx = TRUE` default yields PowerApprox (SOC).  The old assertion
  ## named the bare `Power` class and so pinned the PowCone3D defect; because
  ## PowerApprox extends Power it could never have caught it.
  expect_s3_class(s, "CVXR::PowerApprox")
  expect_s3_class(s, "CVXR::Power")
  expect_equal(as.numeric(s@p_used), 0.5)
  expect_true(is_concave(s))
})

## ── Math group: cumsum ───────────────────────────────────────────

## @cvxpy NONE
test_that("cumsum(x) creates Cumsum atom via Math dispatch", {
  x <- Variable(3)
  cs <- cumsum(x)
  expect_s3_class(cs, "CVXR::Cumsum")
  expect_equal(cs@shape, c(3L, 1L))
})

## ── Math group: sign(x) routes to Sign atom ──────────────────────

## @cvxpy NONE
## sign(x) used to abort with "use expr_sign" -- that error has been
## replaced by routing through the Math handler to the Sign DQCP atom
## (v1.8.2-1, sign-atom branch).  expr_sign() is still the
## curvature-query generic; they are distinct concepts.
test_that("sign(x) dispatches to the Sign atom for CVXR expressions", {
  x <- Variable(3)
  s <- sign(x)
  expect_s3_class(s, "CVXR::Sign")
  expect_equal(s@shape, c(3L, 1L))
})

## ── Math group: log2, log10, log1p ───────────────────────────────

## @cvxpy NONE
test_that("log2(x) works via Math dispatch", {
  x <- Variable(3)
  l2 <- log2(x)
  ## Should be Log(x) / Constant(log(2))
  expect_s3_class(l2, "CVXR::DivExpression")
})

## @cvxpy NONE
test_that("log10(x) works via Math dispatch", {
  x <- Variable(3)
  l10 <- log10(x)
  expect_s3_class(l10, "CVXR::DivExpression")
})

## @cvxpy NONE
test_that("log1p(x) works via Math dispatch", {
  x <- Variable(3)
  l1p <- log1p(x)
  ## Should be Log(x + 1)
  expect_s3_class(l1p, "CVXR::Log")
})

## ── Math group: unsupported ──────────────────────────────────────

## @cvxpy NONE
test_that("Unsupported Math functions give clear error", {
  ## NOTE: sin/cos/tan (and sinh/tanh/asinh/atanh) became supported DNLP atoms
  ## in CVXPY 1.9.0 (see test-dnlp-atoms.R).  These remain unsupported: cosh has
  ## no atom upstream, and the inverse-trig / gamma functions are not modeled.
  x <- Variable(3)
  expect_error(cosh(x), "not supported")
  expect_error(acos(x), "not supported")
  expect_error(gamma(x), "not supported")
})

## @cvxpy NONE
test_that("cummax(x) creates Cummax atom via Math dispatch", {
  x <- Variable(3)
  cm <- cummax(x)
  expect_s3_class(cm, "CVXR::Cummax")
  expect_equal(cm@shape, c(3L, 1L))
  expect_true(is_convex(cm))
})

## @cvxpy NONE
test_that("cumprod(x) creates Cumprod atom via Math dispatch", {
  x <- Variable(3)
  cp <- cumprod(x)
  expect_s3_class(cp, "CVXR::Cumprod")
  expect_equal(cp@shape, c(3L, 1L))
})

## ── Summary group: sum ───────────────────────────────────────────

## @cvxpy NONE
test_that("sum(x) creates SumEntries atom via Summary dispatch", {
  x <- Variable(c(3, 2))
  s <- sum(x)
  expect_s3_class(s, "CVXR::SumEntries")
  expect_equal(s@shape, c(1L, 1L))
  expect_true(is_convex(s))
  expect_true(is_concave(s))
})

## @cvxpy NONE
test_that("sum(x, y) adds then sums", {
  x <- Variable(3)
  y <- Variable(3)
  s <- sum(x, y)
  ## Should be SumEntries of (x + y)
  expect_s3_class(s, "CVXR::SumEntries")
})

## ── Summary group: max ───────────────────────────────────────────

## @cvxpy NONE
test_that("max(x) creates MaxEntries atom via Summary dispatch", {
  x <- Variable(3)
  m <- max(x)
  expect_s3_class(m, "CVXR::MaxEntries")
  expect_equal(m@shape, c(1L, 1L))
  expect_true(is_convex(m))
})

## @cvxpy NONE
test_that("max(x, y) creates Maximum (elementwise) via Summary dispatch", {
  x <- Variable(3)
  y <- Variable(3)
  m <- max(x, y)
  ## With multiple args, max should route to Maximum
  expect_s3_class(m, "CVXR::Maximum")
})

## ── Summary group: min ───────────────────────────────────────────

## @cvxpy NONE
test_that("min(x) creates MinEntries atom via Summary dispatch", {
  x <- Variable(3)
  m <- min(x)
  expect_s3_class(m, "CVXR::MinEntries")
  expect_equal(m@shape, c(1L, 1L))
  expect_true(is_concave(m))
})

## @cvxpy NONE
test_that("min(x, y) creates Minimum (elementwise) via Summary dispatch", {
  x <- Variable(3)
  y <- Variable(3)
  m <- min(x, y)
  expect_s3_class(m, "CVXR::Minimum")
})

## ── Summary group: unsupported ───────────────────────────────────

## @cvxpy NONE
test_that("prod(x) creates Prod atom via Summary dispatch", {
  x <- Variable(3)
  p <- prod(x)
  expect_s3_class(p, "CVXR::Prod")
  expect_equal(p@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Unsupported Summary functions give clear error", {
  x <- Variable(3)
  expect_error(range(x), "not supported")
})

## ── Ops: ^ ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("x^p creates Power atom via Ops dispatch", {
  x <- Variable(3)
  p <- x^2
  expect_s3_class(p, "CVXR::Power")
  expect_equal(p@p_used, 2)
  expect_true(is_convex(p))
})

## @cvxpy NONE
test_that("x^0.5 is concave via Ops dispatch", {
  x <- Variable(3)
  p <- x^0.5
  expect_s3_class(p, "CVXR::Power")
  expect_true(is_concave(p))
})

## ── Nested compositions via natural R syntax ─────────────────────

## @cvxpy NONE
test_that("Nested: exp(sum(x)) is convex", {
  x <- Variable(3)
  e <- exp(sum(x))
  expect_true(is_convex(e))
  ## sum(x) is affine, exp is convex+increasing → composition is convex
})

## @cvxpy NONE
test_that("Nested: sum(abs(x)) is convex", {
  x <- Variable(3)
  s <- sum(abs(x))
  expect_true(is_convex(s))
})

## @cvxpy NONE
test_that("Nested: log(sum(exp(x))) is convex — DCP violation", {
  x <- Variable(3)
  ## log is concave, but sum(exp(x)) is convex not concave
  ## so log(sum(exp(x))) violates DCP
  l <- log(sum(exp(x)))
  ## This should NOT be recognized as convex by DCP rules
  ## (even though log-sum-exp IS convex, DCP can't verify it)
  expect_false(is_convex(l))
  expect_false(is_concave(l))
})
