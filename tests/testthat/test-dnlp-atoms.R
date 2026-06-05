# Tests for the DNLP smooth atoms (CVXPY 1.9.0): trig (sin/cos/tan),
# hyperbolic (sinh/tanh/asinh/atanh), and normcdf.
#
# These are `now`-bucket atoms. They are smooth but neither convex nor
# concave, so they are NOT DCP -- they participate only on the disciplined-
# nonlinear-programming (DNLP) path. The dnlp2smooth reduction + NLP solving
# chain have landed, so both atom-level behavior and the end-to-end
# `nlp = TRUE` solves are active here.
#
# Expected values cross-checked against R's own sin/cos/.../pnorm (which match
# numpy/scipy used by CVXPY) and the curvature/sign/monotonicity declared in
# the cited source.
#
# CVXPY SOURCE: atoms/elementwise/trig.py, atoms/elementwise/hyperbolic.py,
#   atoms/elementwise/normcdf.py

library(CVXR)

# ===================================================================
# Construction + curvature/sign/monotonicity/smoothness
# ===================================================================

## @cvxpy NONE
test_that("trig atoms: smooth, not convex/concave, sign unknown, not monotone", {
  x <- Variable(3)
  for (a in list(Sin(x), Cos(x), Tan(x))) {
    expect_equal(a@shape, c(3L, 1L))
    expect_true(is_atom_smooth(a))
    expect_false(is_atom_convex(a))
    expect_false(is_atom_concave(a))
    expect_false(is_incr(a, 1L))
    expect_false(is_decr(a, 1L))
    sgn <- sign_from_args(a)
    expect_false(sgn$is_nonneg)
    expect_false(sgn$is_nonpos)
  }
})

## @cvxpy NONE
test_that("hyperbolic atoms: smooth, not convex/concave, increasing, sign unknown", {
  x <- Variable(3)
  for (a in list(Sinh(x), Tanh(x), Asinh(x), Atanh(x))) {
    expect_equal(a@shape, c(3L, 1L))
    expect_true(is_atom_smooth(a))
    expect_false(is_atom_convex(a))
    expect_false(is_atom_concave(a))
    expect_true(is_incr(a, 1L))
    expect_false(is_decr(a, 1L))
    sgn <- sign_from_args(a)
    expect_false(sgn$is_nonneg)
    expect_false(sgn$is_nonpos)
  }
})

## @cvxpy NONE
test_that("normcdf: smooth, nonneg, increasing, not convex/concave", {
  x <- Variable(3)
  a <- Normcdf(x)
  expect_equal(a@shape, c(3L, 1L))
  expect_true(is_atom_smooth(a))
  expect_false(is_atom_convex(a))
  expect_false(is_atom_concave(a))
  expect_true(is_incr(a, 1L))
  expect_false(is_decr(a, 1L))
  sgn <- sign_from_args(a)
  expect_true(sgn$is_nonneg)
  expect_false(sgn$is_nonpos)
})

# ===================================================================
# DNLP / DCP classification (ties into the Step-2 predicate cluster)
# ===================================================================

## @cvxpy NONE
test_that("smooth atoms are DNLP but not DCP", {
  x <- Variable(2)
  for (e in list(Sin(x), Cos(x), Tan(x), Sinh(x), Tanh(x),
                 Asinh(x), Atanh(x), Normcdf(x))) {
    ## smooth(arg=leaf) => linearizable both ways => is_smooth => is_dnlp
    expect_true(is_smooth(e))
    expect_true(is_dnlp(e))
    ## neither convex nor concave => not DCP
    expect_false(is_dcp(e))
  }
})

# ===================================================================
# Numeric evaluation (matches numpy/scipy via R's base functions)
# ===================================================================

## @cvxpy NONE
test_that("trig: numeric evaluation", {
  x <- Variable(2)
  val <- matrix(c(0, 1), ncol = 1)
  expect_equal(numeric_value(Sin(x), list(val)), matrix(sin(c(0, 1)), ncol = 1))
  expect_equal(numeric_value(Cos(x), list(val)), matrix(cos(c(0, 1)), ncol = 1))
  expect_equal(numeric_value(Tan(x), list(val)), matrix(tan(c(0, 1)), ncol = 1))
})

## @cvxpy NONE
test_that("hyperbolic: numeric evaluation", {
  x <- Variable(2)
  val <- matrix(c(0, 1), ncol = 1)
  vatanh <- matrix(c(0, 0.5), ncol = 1)  # atanh domain is (-1, 1)
  expect_equal(numeric_value(Sinh(x), list(val)),  matrix(sinh(c(0, 1)), ncol = 1))
  expect_equal(numeric_value(Tanh(x), list(val)),  matrix(tanh(c(0, 1)), ncol = 1))
  expect_equal(numeric_value(Asinh(x), list(val)), matrix(asinh(c(0, 1)), ncol = 1))
  expect_equal(numeric_value(Atanh(x), list(vatanh)), matrix(atanh(c(0, 0.5)), ncol = 1))
})

## @cvxpy NONE
test_that("normcdf: numeric evaluation equals pnorm", {
  x <- Variable(2)
  val <- matrix(c(0, 1), ncol = 1)
  expect_equal(numeric_value(Normcdf(x), list(val)), matrix(pnorm(c(0, 1)), ncol = 1))
})

# ===================================================================
# Gradients (implemented for sin/cos/tan/normcdf; NotImplemented for hyperbolic)
# ===================================================================

## @cvxpy NONE
test_that("trig gradients: diagonal of derivative at the point", {
  x <- Variable(2)
  val <- matrix(c(0, 1), ncol = 1)
  expect_equal(diag(as.matrix(.grad(Sin(x), list(val))[[1L]])), cos(c(0, 1)))
  expect_equal(diag(as.matrix(.grad(Cos(x), list(val))[[1L]])), -sin(c(0, 1)))
  expect_equal(diag(as.matrix(.grad(Tan(x), list(val))[[1L]])), 1 / cos(c(0, 1))^2)
})

## @cvxpy NONE
test_that("normcdf gradient equals dnorm on the diagonal", {
  x <- Variable(2)
  val <- matrix(c(0, 1), ncol = 1)
  expect_equal(diag(as.matrix(.grad(Normcdf(x), list(val))[[1L]])), dnorm(c(0, 1)))
})

## @cvxpy NONE
test_that("hyperbolic gradients are not implemented (mirrors CVXPY)", {
  x <- Variable(2)
  val <- matrix(c(0, 0.5), ncol = 1)
  expect_error(.grad(Sinh(x), list(val)),  "not implemented")
  expect_error(.grad(Tanh(x), list(val)),  "not implemented")
  expect_error(.grad(Asinh(x), list(val)), "not implemented")
  expect_error(.grad(Atanh(x), list(val)), "not implemented")
})

# ===================================================================
# Domains
# ===================================================================

## @cvxpy NONE
test_that("tan domain is the two box constraints -pi/2 <= x <= pi/2", {
  x <- Variable(2)
  d <- atom_domain(Tan(x))
  expect_equal(length(d), 2L)
  expect_true(all(vapply(d, function(c) S7_inherits(c, Inequality), logical(1))))
})

## @cvxpy NONE
test_that("atanh domain aborts: open interval needs strict inequalities", {
  x <- Variable(2)
  expect_error(atom_domain(Atanh(x)), "strict inequalities")
})

# ===================================================================
# Standard-R Math-group dispatch routes to the atoms
# ===================================================================

## @cvxpy NONE
test_that("base R trig/hyperbolic functions build the corresponding atoms", {
  x <- Variable(3)
  expect_true(S7_inherits(sin(x),   Sin))
  expect_true(S7_inherits(cos(x),   Cos))
  expect_true(S7_inherits(tan(x),   Tan))
  expect_true(S7_inherits(sinh(x),  Sinh))
  expect_true(S7_inherits(tanh(x),  Tanh))
  expect_true(S7_inherits(asinh(x), Asinh))
  expect_true(S7_inherits(atanh(x), Atanh))
  expect_true(S7_inherits(normcdf(x), Normcdf))
})

## @cvxpy NONE
test_that("cosh has no atom (CVXPY 1.9.0 defines none) and is rejected", {
  x <- Variable(3)
  expect_error(cosh(x), "not supported")
})

# ===================================================================
# End-to-end nlp solves -- pending the DNLP solving chain (Wave-1 steps 4-7)
# ===================================================================

## @cvxpy nlp_tests/test_hyperbolic.py::TestHyperbolic::test_sinh
test_that("nlp solve of a sinh objective (DNLP path)", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## sinh is smooth and strictly increasing, so min sinh(x) s.t. x >= 0.5 is
  ## attained at the lower bound: x* = 0.5, f* = sinh(0.5). Genuinely DNLP
  ## (sinh is neither convex nor concave), exercising the smooth-atom solve path.
  x <- Variable()
  value(x) <- 1
  prob <- Problem(Minimize(sinh(x)), list(x >= 0.5))
  val <- psolve(prob, solver = "UNO", nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, sinh(0.5), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_hyperbolic.py::TestHyperbolic::test_atanh
test_that("nlp solve of an atanh objective (DNLP path)", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  ## Mirrors CVXPY test_atanh: atanh's domain is (-1, 1), kept satisfied by
  ## composing it with logistic(x * 0.1) (whose values stay well inside the
  ## domain at the iterates) rather than by any bound -- the solver simply
  ## evaluates the oracle. min sum(atanh(logistic(x * 0.1))) s.t. x >= 0.1,
  ## sum(x) == 10.
  n <- 10
  x <- Variable(n)
  prob <- Problem(Minimize(sum_entries(atanh(logistic(x * 0.1)))),
                  list(x >= 0.1, sum_entries(x) == 10))
  val <- psolve(prob, solver = "UNO", nlp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 9.602686, tolerance = 1e-3)
})
