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
  ## sqrt.py:24 is `power(x, Fraction(1, 2))` -- the WRAPPER, whose
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

## numpy's log has no base; this is base R's log(x, base).
## @cvxpy NONE
test_that("log(x, base) honors base", {
  v <- matrix(c(2, 3, 4.5, 5, 1, 7), 2, 3)
  expect_equal(value(log(Constant(v), 2)), log(v, 2))
  expect_equal(value(log(Constant(v), base = 10)), log(v, base = 10))
  expect_equal(value(log(Constant(v))), log(v))
  x <- Variable(3)
  expect_true(is_concave(log(x, 2)))
  ## a base below 1 is a negative divisor and flips the curvature
  expect_true(is_convex(log(x, 0.5)))
  expect_error(log(x, 1), "base")
  expect_error(log(x, 0), "base")
  expect_error(log(x, -2), "base")
  expect_error(log(x, c(2, 3)), "base")
  expect_error(log(x, Variable()), "base")
})

## ── Math group: unsupported ──────────────────────────────────────

## elementwise/log1p.py:25-41: log1p is its own atom, evaluated with
## scipy.special.log1p, with the sign of its argument.
## @cvxpy NONE
test_that("log1p(x) is the Log1p atom", {
  x <- Variable(3, nonneg = TRUE)
  expect_true(S7::S7_inherits(log1p(x), Log1p))
  expect_true(is_nonneg(log1p(x)))
  expect_true(is_concave(log1p(x)))
  v <- c(1e-10, 1e-3, 2)
  expect_equal(as.numeric(value(log1p(Constant(v)))), log1p(v))
  ## log(1 + 1e-10) is off by a relative 8e-8, which expect_equal() cannot
  ## see at this magnitude (it compares absolutely); Log1p calls log1p()
  expect_identical(as.numeric(value(log1p(Constant(1e-10)))), log1p(1e-10))
  prob <- Problem(Maximize(sum(log1p(x))), list(sum(x) <= 1))
  expect_equal(psolve(prob, solver = "CLARABEL"), 3 * log1p(1 / 3),
               tolerance = 1e-6)
})

## stats.py:35 mean(x, axis=None, keepdims=False). mean() is a generic, so
## trim/na.rm are mean.default's; on an expression they are errors and a
## positional second argument is axis.
## @cvxpy NONE
test_that("mean(expr) is cvxr_mean; trim and na.rm are errors", {
  m <- matrix(c(2, 3, 4.5, 5, 1, 7), 2, 3)
  X <- Constant(m)
  expect_equal(as.numeric(value(mean(X))), mean(m))
  expect_equal(as.numeric(value(mean(X, 1))), rowMeans(m))
  expect_equal(as.numeric(value(mean(X, axis = 2))), colMeans(m))
  x <- Variable(c(2, 3))
  expect_error(mean(x, trim = 0.1), "cvxr_mean")
  expect_error(mean(x, na.rm = TRUE), "cvxr_mean")
  expect_error(mean(x, 0.1), "axis")
  expect_error(mean(x, 1.5), "axis")
  expect_error(mean(x, c(1, 2)), "axis")
  expect_error(mean(x, x), "axis")
  ## plain objects are base R's mean
  expect_identical(mean(c(1, 2, 30), trim = 0.4), 2)
})

## power.py:183 float(self.p.value): a constant exponent must be a scalar.
## @cvxpy NONE
test_that("A non-scalar constant exponent is an error", {
  x <- Variable(2)
  expect_error(x^c(1, 2), "scalar")
  expect_error(x^numeric(0), "scalar")
  expect_error(power(x, c(1, 2)), "scalar")
  expect_error(power(x, c(1, 2), approx = FALSE), "scalar")
  expect_error(power(x, Constant(c(1, 2))), "scalar")
  expect_error(2^Constant(c(1, 2)), "scalar")
  expect_error(Constant(2)^c(1, 2), "scalar")
  ## a 1 x 1 exponent, a variable exponent and a Parameter are unchanged
  expect_true(S7_inherits(x^matrix(2), Power))
  expect_true(S7_inherits(power(x, Constant(matrix(3)), approx = FALSE), Power))
  expect_false(S7_inherits(2^x, Power))
  expect_true(S7_inherits(x^Parameter(), Power))
})

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
test_that("sum(x, y, ...) is the sum of the entry sums, as base::sum", {
  ## base::sum(x, y) adds every entry of every operand; shapes need not
  ## agree. sum(x, y) was sum(x + y) with broadcasting, so sum(x, 2) for a
  ## 2 x 3 x added 2 six times.
  x <- Variable(c(2L, 3L))
  y <- Variable(c(3L, 1L))
  M <- matrix(c(1, -4, 3, 2, -5, 6), 2L, 3L)
  N <- matrix(c(10, 20, 30), 3L, 1L)
  value(x) <- M
  value(y) <- N
  s <- sum(x, y)
  expect_equal(s@shape, c(1L, 1L))
  expect_equal(as.numeric(value(s)), sum(M, N))
  expect_equal(as.numeric(value(sum(x, 2))), sum(M, 2))
  expect_equal(as.numeric(value(sum(x, y, 2))), sum(M, N, 2))
  expect_true(is_affine(s))
  ## prod likewise: the product of the entry products, not prod(x * y).
  expect_equal(as.numeric(value(prod(x, 2))), prod(M, 2))
  expect_equal(as.numeric(value(prod(x, y))), prod(M, N))
  expect_equal(prod(x, y)@shape, c(1L, 1L))
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
test_that("max(x, y, ...) is one scalar over all operands, as base::max", {
  ## base::max(x, y) is the largest entry of all operands, a length-1 result;
  ## the elementwise maximum is pmax(), max_elemwise() here. max(x, y) used
  ## to be the elementwise Maximum(x, y).
  x <- Variable(c(2L, 3L))
  y <- Variable(c(3L, 1L))
  M <- matrix(c(1, -4, 3, 2, -5, 6), 2L, 3L)
  N <- matrix(c(10, -20, 3), 3L, 1L)
  value(x) <- M
  value(y) <- N
  m <- max(x, y)
  expect_s3_class(m, "CVXR::Maximum")
  expect_equal(m@shape, c(1L, 1L))
  expect_equal(as.numeric(value(m)), max(M, N))
  expect_true(is_convex(m))
  expect_equal(max(x, 0)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(max(x, 0))), max(M, 0))
  expect_equal(as.numeric(value(max(x, y, 100))), 100)
  ## The elementwise maximum is still max_elemwise().
  expect_equal(max_elemwise(x, 0)@shape, c(2L, 3L))
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
test_that("min(x, y, ...) is one scalar over all operands, as base::min", {
  x <- Variable(c(2L, 3L))
  y <- Variable(c(3L, 1L))
  M <- matrix(c(1, -4, 3, 2, -5, 6), 2L, 3L)
  N <- matrix(c(10, -20, 3), 3L, 1L)
  value(x) <- M
  value(y) <- N
  m <- min(x, y)
  expect_s3_class(m, "CVXR::Minimum")
  expect_equal(m@shape, c(1L, 1L))
  expect_equal(as.numeric(value(m)), min(M, N))
  expect_true(is_concave(m))
  expect_equal(min(x, 0)@shape, c(1L, 1L))
  expect_equal(as.numeric(value(min(x, 0))), min(M, 0))
  expect_equal(min_elemwise(x, 0)@shape, c(2L, 3L))
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
test_that("sum/max/min/prod honor axis and keepdims as CVXPY's do", {
  ## A named axis was taken as one more operand: sum(x, axis = 2L) was
  ## sum(x + 2) and max(x, axis = 2L) was the elementwise max with 2.
  x <- Variable(c(2L, 3L))
  M <- matrix(c(1, -4, 3, 2, -5, 6), 2L, 3L)
  value(x) <- M
  cases <- list(sum = list(sum, sum_entries, base::sum),
                max = list(max, max_entries, base::max),
                min = list(min, min_entries, base::min),
                prod = list(prod, function(e, ...) Prod(e, ...), base::prod))
  for (nm in names(cases)) {
    f <- cases[[nm]][[1L]]; g <- cases[[nm]][[2L]]; b <- cases[[nm]][[3L]]
    for (ax in 1:2) {
      e <- f(x, axis = ax)
      expect_equal(e@shape, g(x, axis = ax)@shape, info = paste(nm, ax))
      expect_equal(as.numeric(value(e)), as.numeric(apply(M, ax, b)),
                   info = paste(nm, ax))
      ek <- f(x, axis = ax, keepdims = TRUE)
      expect_equal(ek@shape, g(x, axis = ax, keepdims = TRUE)@shape,
                   info = paste(nm, ax, "keepdims"))
    }
    expect_equal(as.numeric(value(f(x, keepdims = TRUE))), b(M), info = nm)
  }
})

## @cvxpy NONE
test_that("sum/max/min/prod take named operands as base R does; axis needs one operand", {
  ## base::sum(..., na.rm) labels a `...` entry by its name and still adds it:
  ## sum(a = 1:3) is 6 and sum(1:3, axes = 2) is 8. Only `axis` and `keepdims`
  ## are options here (the CVXPY keywords); every other entry is an operand.
  x <- Variable(c(2L, 3L))
  y <- Variable(c(2L, 3L))
  M <- matrix(c(1, -4, 3, 2, -5, 6), 2L, 3L)
  value(x) <- M
  value(y) <- 2 * M
  expect_s3_class(sum(a = x), "CVXR::SumEntries")
  expect_equal(as.numeric(value(sum(a = x))), sum(a = M))
  expect_equal(as.numeric(value(do.call(sum, list(a = x, b = y)))), sum(M, 2 * M))
  ## A misspelled option is an operand, as in base R; it means what the
  ## unnamed call means.
  expect_equal(value(sum(x, axes = 2L)), value(sum(x, 2L)))
  expect_s3_class(max(a = x, b = y), "CVXR::Maximum")
  expect_equal(as.numeric(value(max(a = x, b = y))), max(a = M, b = 2 * M))
  expect_equal(as.numeric(value(min(a = x, b = y))), min(a = M, b = 2 * M))
  expect_equal(as.numeric(value(prod(a = x))), prod(a = M))
  ## Operand names do not reach the atom.
  expect_null(names(max(a = x, b = y)@args))
  for (f in list(sum, max, min, prod)) {
    expect_error(f(x, y, axis = 2L), "single expression")
    expect_error(f(a = x, b = y, keepdims = TRUE), "single expression")
  }
  ## Unnamed operands keep the base R meaning.
  expect_equal(sum(x, y)@shape, c(1L, 1L))
  expect_equal(max(x, y)@shape, c(1L, 1L))
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
