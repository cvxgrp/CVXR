# Tests for the Dnlp2Smooth reduction + canonicalizers (CVXPY 1.9.0, Step 4a).
#
# Dnlp2Smooth reduces a disciplined nonlinear program to an equivalent SMOOTH
# program: smooth atoms (exp/log/sin/...) are kept (the diff engine handles
# their chain rule), nonsmooth-but-representable atoms (quad_over_lin, pnorm,
# power, rel_entr, ...) get smooth epigraph reformulations with aux variables
# and equality constraints. The NLP solving chain + Uno backend have landed, so
# the end-to-end nlp=TRUE solves are active (see the test-nlp-*-v19.R family).
#
# CVXPY SOURCE: reductions/dnlp2smooth/dnlp2smooth.py + canonicalizers/*.py

library(CVXR)

# ===================================================================
# Reduction-level behavior
# ===================================================================

## @cvxpy NONE
test_that("Dnlp2Smooth accepts any problem", {
  x <- Variable(2)
  expect_true(reduction_accepts(Dnlp2Smooth(), Problem(Minimize(sum_entries(x)), list())))
})

## @cvxpy NONE
test_that("a smooth chain-rule atom (exp) is kept, no aux constraints", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list())
  res <- reduction_apply(Dnlp2Smooth(), prob)
  newprob <- res[[1L]]
  expect_true(S7_inherits(newprob, Problem))
  expect_equal(length(newprob@constraints), 0L)
})

## @cvxpy NONE
test_that("a Step-3 smooth trig atom (sin) survives the reduction unchanged", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Sin(x))), list())
  expect_true(is_dnlp(prob))
  res <- reduction_apply(Dnlp2Smooth(), prob)
  expect_equal(length(res[[1L]]@constraints), 0L)
})

## @cvxpy NONE
test_that("a nonnegative-domain atom (log) is lifted with one equality", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Log(x))), list())
  res <- reduction_apply(Dnlp2Smooth(), prob)
  ## log_canon introduces a nonneg aux variable t with t == x
  expect_equal(length(res[[1L]]@constraints), 1L)
})

## @cvxpy NONE
test_that("constraint IDs are recorded in the inverse data cons_id_map", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(x >= 1))
  res <- reduction_apply(Dnlp2Smooth(), prob)
  inv <- res[[2L]]
  expect_true(length(ls(inv@cons_id_map)) >= 1L)
})

# ===================================================================
# Canonicalizer unit tests
# ===================================================================

## @cvxpy NONE
test_that("power_canon keeps integer powers, lifts fractional powers", {
  x <- Variable(2)
  p2 <- Power(x, 2)
  r2 <- .smooth_power_canon(p2, p2@args)
  expect_true(S7_inherits(r2[[1L]], Power))
  expect_equal(length(r2[[2L]]), 0L)

  ph <- Power(x, 0.5)
  rh <- .smooth_power_canon(ph, ph@args)
  expect_equal(length(rh[[2L]]), 1L)   # t == x (nonneg aux)
})

## @cvxpy NONE
test_that("log_canon introduces a nonneg aux variable + equality", {
  x <- Variable(2)
  l <- Log(x)
  r <- .smooth_log_canon(l, l@args)
  expect_true(S7_inherits(r[[1L]], Log))
  expect_equal(length(r[[2L]]), 1L)
  expect_true(S7_inherits(r[[1L]]@args[[1L]], Variable))
})

## @cvxpy NONE
test_that("div_canon rewrites division as multiply by reciprocal", {
  x <- Variable(2)
  y <- Variable(2, nonneg = TRUE)
  d <- x / y
  r <- .smooth_div_canon(d, d@args)
  ## one equality (z == denominator); result is an expression
  expect_equal(length(r[[2L]]), 1L)
  expect_true(S7_inherits(r[[1L]], Expression))
})

## @cvxpy NONE
test_that("quad_over_lin canon with constant denominator reuses power_canon", {
  x <- Variable(2)
  qol <- QuadOverLin(x, 2)
  r <- .smooth_quad_over_lin_canon(qol, qol@args)
  ## constant-denominator branch: (1/2)*sum(x^2); p=2 is integer -> no constraints
  expect_equal(length(r[[2L]]), 0L)
  expect_true(S7_inherits(r[[1L]], Expression))
})

## @cvxpy NONE
test_that("log_sum_exp canon introduces t, nonpos v, and two constraints", {
  x <- Variable(3)
  lse <- LogSumExp(x)
  r <- .smooth_log_sum_exp_canon(lse, lse@args)
  expect_true(S7_inherits(r[[1L]], Variable))   # returns t
  expect_equal(length(r[[2L]]), 2L)             # sum(exp(v))==1, v==x-t
})

# ===================================================================
# Step 4b: nonsmooth PWL atoms reuse the dcp2cone canonicalizers
# ===================================================================

## @cvxpy NONE
test_that("PWL atoms are registered for smooth_canonicalize (not the NULL default)", {
  x <- Variable(3)
  expect_false(is.null(smooth_canonicalize(Abs(x), list(x))))
  expect_false(is.null(smooth_canonicalize(MaxEntries(x), list(x))))
  expect_false(is.null(smooth_canonicalize(MinEntries(x), list(x))))
  expect_false(is.null(smooth_canonicalize(Norm1(x), list(x))))
  expect_false(is.null(smooth_canonicalize(NormInf(x), list(x))))
  expect_false(is.null(smooth_canonicalize(SumLargest(x, 2), list(x))))
  m <- Maximum(x, 0)
  expect_false(is.null(smooth_canonicalize(m, m@args)))
  mn <- Minimum(x, 0)
  expect_false(is.null(smooth_canonicalize(mn, mn@args)))
})

## @cvxpy NONE
test_that("smooth_canonicalize on abs matches the dcp2cone abs_canon epigraph", {
  x <- Variable(3)
  a <- Abs(x)
  r <- smooth_canonicalize(a, list(x))
  ## abs epigraph: t with t >= x, t >= -x  (two constraints)
  expect_equal(length(r[[2L]]), length(abs_canon(a, list(x))[[2L]]))
})

## @cvxpy NONE
test_that("Dnlp2Smooth lifts a nonsmooth abs() objective to an epigraph", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(Abs(x))), list())
  res <- reduction_apply(Dnlp2Smooth(), prob)
  expect_true(length(res[[1L]]@constraints) >= 1L)
})

# ===================================================================
# End-to-end nlp solve -- pending the DNLP solving chain (steps 5-7)
# ===================================================================

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_simple_smooth
test_that("nlp=TRUE solve of a DNLP problem", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## A DCP (hence DNLP) problem routed through Dnlp2Smooth + Uno.
  ## min sum(exp(x))  s.t.  sum(x) == 1.5  =>  x_i = 0.5, f* = 3 exp(0.5).
  x <- Variable(3)
  value(x) <- rep(0, 3)
  prob <- Problem(Minimize(sum_entries(Exp(x))), list(sum_entries(x) == 1.5))
  val <- psolve(prob, solver = "UNO", nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 3 * exp(0.5), tolerance = 1e-4)
})
