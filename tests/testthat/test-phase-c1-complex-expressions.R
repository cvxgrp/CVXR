## Tests for Phase C1: Complex expression-level infrastructure
## Conj_, Real_, Imag_ atoms, Wrap classes, Complex group handler, expr_H

# ── Conj_ atom ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Conj_ creates correct shape", {
  x <- Variable(c(3L, 2L))
  cx <- Conj_(x)
  expect_equal(cx@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Conj_ numeric_value uses R's Conj()", {
  v <- matrix(c(1+2i, 3+4i), 2, 1)
  c1 <- Constant(v)
  cx <- Conj_(c1)
  result <- numeric_value(cx, list(v))
  expect_equal(result, Conj(v))
})

## @cvxpy NONE
test_that("Conj_ numeric_value is no-op for real", {
  v <- matrix(c(1, 2, 3), 3, 1)
  c1 <- Constant(v)
  cx <- Conj_(c1)
  result <- numeric_value(cx, list(v))
  expect_equal(result, v)
})

## @cvxpy NONE
test_that("Conj_ is not monotone", {
  x <- Variable(c(2L, 2L))
  cx <- Conj_(x)
  expect_false(is_incr(cx, 1L))
  expect_false(is_decr(cx, 1L))
})

## @cvxpy NONE
test_that("Conj_ delegates is_symmetric to arg", {
  x <- Variable(c(3L, 3L), symmetric = TRUE)
  cx <- Conj_(x)
  expect_true(is_symmetric(cx))

  y <- Variable(c(3L, 2L))
  cy <- Conj_(y)
  expect_false(is_symmetric(cy))
})

## @cvxpy NONE
test_that("Conj_ delegates is_hermitian to arg", {
  x <- Variable(c(3L, 3L), hermitian = TRUE)
  cx <- Conj_(x)
  expect_true(is_hermitian(cx))

  y <- Variable(c(3L, 3L))
  cy <- Conj_(y)
  expect_false(is_hermitian(cy))
})

## @cvxpy NONE
test_that("Conj_ graph_implementation is identity passthrough", {
  x <- Variable(c(2L, 1L))
  cx <- Conj_(x)
  dummy_obj <- list("fake_linop")
  result <- graph_implementation(cx, list(dummy_obj), c(2L, 1L))
  expect_identical(result[[1L]], dummy_obj)
  expect_length(result[[2L]], 0L)
})

## @cvxpy NONE
test_that("Conj_ is affine", {
  x <- Variable(c(2L, 1L))
  cx <- Conj_(x)
  expect_true(is_convex(cx))
  expect_true(is_concave(cx))
  expect_true(is_affine(cx))
})

# ── Real_ atom ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Real_ creates correct shape", {
  x <- Variable(c(3L, 2L))
  rx <- Real_(x)
  expect_equal(rx@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Real_ numeric_value uses R's Re()", {
  v <- matrix(c(1+2i, 3+4i), 2, 1)
  c1 <- Constant(v)
  rx <- Real_(c1)
  result <- numeric_value(rx, list(v))
  expect_equal(result, Re(v))
})

## @cvxpy NONE
test_that("Real_ is never complex or imaginary", {
  x <- Variable(c(2L, 1L), complex = TRUE)
  rx <- Real_(x)
  expect_false(is_complex(rx))
  expect_false(is_imag(rx))
})

## @cvxpy NONE
test_that("Real_ is_symmetric delegates to arg is_hermitian", {
  h <- Variable(c(3L, 3L), hermitian = TRUE)
  rx <- Real_(h)
  expect_true(is_symmetric(rx))

  y <- Variable(c(3L, 3L))
  ry <- Real_(y)
  expect_false(is_symmetric(ry))
})

## @cvxpy NONE
test_that("Real_ is affine", {
  x <- Variable(c(2L, 1L))
  rx <- Real_(x)
  expect_true(is_affine(rx))
})

# ── Imag_ atom ───────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Imag_ creates correct shape", {
  x <- Variable(c(3L, 2L))
  ix <- Imag_(x)
  expect_equal(ix@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Imag_ numeric_value uses R's Im()", {
  v <- matrix(c(1+2i, 3+4i), 2, 1)
  c1 <- Constant(v)
  ix <- Imag_(c1)
  result <- numeric_value(ix, list(v))
  expect_equal(result, Im(v))
})

## @cvxpy NONE
test_that("Imag_ is never complex or imaginary", {
  x <- Variable(c(2L, 1L), complex = TRUE)
  ix <- Imag_(x)
  expect_false(is_complex(ix))
  expect_false(is_imag(ix))
})

## @cvxpy NONE
test_that("Imag_ is_symmetric: Im(Hermitian) is skew-symmetric, not symmetric", {
  ## CVXPY v1.8.2 fix: is_symmetric now delegates to arg is_symmetric (not is_hermitian).
  ## Im(Hermitian) is skew-symmetric. Im(symmetric real) is zero (symmetric).
  h <- Variable(c(3L, 3L), hermitian = TRUE)
  ix <- Imag_(h)
  expect_false(is_symmetric(ix))  # skew-symmetric, not symmetric

  s <- Variable(c(3L, 3L), symmetric = TRUE)
  expect_true(is_symmetric(Imag_(s)))  # Im of real symmetric is zero (symmetric)
})

## @cvxpy NONE
test_that("Imag_ is affine", {
  x <- Variable(c(2L, 1L))
  ix <- Imag_(x)
  expect_true(is_affine(ix))
})

# ── Complex S3 group handler: Re(), Im(), Conj() ────────────────────

## @cvxpy NONE
test_that("Re() dispatches to Real_ for CVXR expressions", {
  x <- Variable(c(2L, 1L))
  rx <- Re(x)
  expect_true(S7_inherits(rx, Real_))
  expect_equal(rx@shape, c(2L, 1L))
})

## @cvxpy NONE
test_that("Im() dispatches to Imag_ for CVXR expressions", {
  x <- Variable(c(2L, 1L))
  ix <- Im(x)
  expect_true(S7_inherits(ix, Imag_))
})

## @cvxpy NONE
test_that("Conj() dispatches to Conj_ for CVXR expressions", {
  x <- Variable(c(2L, 1L))
  cx <- Conj(x)
  expect_true(S7_inherits(cx, Conj_))
})

## @cvxpy NONE
test_that("Mod() dispatches to Abs for CVXR expressions", {
  x <- Variable(c(2L, 1L))
  mx <- Mod(x)
  expect_true(S7_inherits(mx, Abs))
})

## @cvxpy NONE
test_that("Arg() errors for CVXR expressions", {
  x <- Variable(c(2L, 1L))
  expect_error(Arg(x), "Arg.*not supported")
})

# ── Complex handler with Constants ───────────────────────────────────

## @cvxpy NONE
test_that("Re() on complex Constant creates Real_ atom", {
  c1 <- Constant(1+2i)
  rc <- Re(c1)
  expect_true(S7_inherits(rc, Real_))
  expect_false(is_complex(rc))
})

## @cvxpy NONE
test_that("Im() on complex Constant creates Imag_ atom", {
  c1 <- Constant(1+2i)
  ic <- Im(c1)
  expect_true(S7_inherits(ic, Imag_))
  expect_false(is_complex(ic))
})

## @cvxpy NONE
test_that("Conj() on real Constant creates Conj_ (no-op)", {
  c1 <- Constant(5)
  cc <- Conj(c1)
  expect_true(S7_inherits(cc, Conj_))
})

# ── expr_H (conjugate-transpose) ────────────────────────────────────

## @cvxpy NONE
test_that("expr_H on real expression is just transpose", {
  x <- Variable(c(3L, 2L))
  h <- expr_H(x)
  ## Real → just Transpose
  expect_true(S7_inherits(h, Transpose))
  expect_equal(h@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("expr_H on complex expression is Conj(t(x))", {
  x <- Variable(c(3L, 2L), complex = TRUE)
  h <- expr_H(x)
  ## Complex → Conj_(Transpose(x))
  expect_true(S7_inherits(h, Conj_))
  inner <- h@args[[1L]]
  expect_true(S7_inherits(inner, Transpose))
  expect_equal(h@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("expr_H on hermitian variable equals self", {
  ## H is square and Hermitian → H^H == H
  ## (structurally, expr_H returns Conj_(t(H)) which is a different object,
  ## but has the same shape)
  H <- Variable(c(3L, 3L), hermitian = TRUE)
  hh <- expr_H(H)
  expect_equal(hh@shape, c(3L, 3L))
})

# ── Wrap base class ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("Wrap is identity passthrough", {
  x <- Variable(c(3L, 2L))
  w <- Wrap(x)
  expect_equal(w@shape, c(3L, 2L))
  expect_true(is_affine(w))
})

## @cvxpy NONE
test_that("Wrap numeric_value returns input", {
  v <- matrix(1:6, 3, 2)
  w <- Wrap(Variable(c(3L, 2L)))
  result <- numeric_value(w, list(v))
  expect_equal(result, v)
})

## @cvxpy NONE
test_that("Wrap delegates is_complex to arg", {
  x <- Variable(c(2L, 1L), complex = TRUE)
  w <- Wrap(x)
  expect_true(is_complex(w))

  y <- Variable(c(2L, 1L))
  w2 <- Wrap(y)
  expect_false(is_complex(w2))
})

# ── nonneg_wrap ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("nonneg_wrap asserts nonneg", {
  x <- Variable(c(3L, 1L))
  w <- nonneg_wrap(x)
  expect_true(is_nonneg(w))
  expect_false(is_nonpos(w))
})

# ── nonpos_wrap ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("nonpos_wrap asserts nonpos", {
  x <- Variable(c(3L, 1L))
  w <- nonpos_wrap(x)
  expect_true(is_nonpos(w))
  expect_false(is_nonneg(w))
})

# ── psd_wrap ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("psd_wrap validates square input", {
  x <- Variable(c(3L, 2L))
  expect_error(psd_wrap(x), "square")
})

## @cvxpy NONE
test_that("psd_wrap asserts PSD, not NSD", {
  x <- Variable(c(3L, 3L))
  w <- psd_wrap(x)
  expect_true(is_psd(w))
  expect_false(is_nsd(w))
})

## @cvxpy NONE
test_that("psd_wrap: real arg is symmetric, complex arg is not", {
  x <- Variable(c(3L, 3L))
  w <- psd_wrap(x)
  expect_true(is_symmetric(w))
  expect_true(is_hermitian(w))

  ## Complex input: Hermitian but NOT symmetric (CVXPY convention)
  z <- Variable(c(3L, 3L), complex = TRUE)
  w2 <- psd_wrap(z)
  expect_false(is_symmetric(w2))
  expect_true(is_hermitian(w2))
})

# ── symmetric_wrap ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("symmetric_wrap validates real square", {
  x <- Variable(c(3L, 2L))
  expect_error(symmetric_wrap(x), "square")

  z <- Variable(c(3L, 3L), complex = TRUE)
  expect_error(symmetric_wrap(z), "real")
})

## @cvxpy NONE
test_that("symmetric_wrap asserts symmetric and Hermitian", {
  x <- Variable(c(3L, 3L))
  w <- symmetric_wrap(x)
  expect_true(is_symmetric(w))
  expect_true(is_hermitian(w))
})

# ── hermitian_wrap ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("hermitian_wrap validates square", {
  x <- Variable(c(3L, 2L))
  expect_error(hermitian_wrap(x), "square")
})

## @cvxpy NONE
test_that("hermitian_wrap asserts Hermitian", {
  x <- Variable(c(3L, 3L))
  w <- hermitian_wrap(x)
  expect_true(is_hermitian(w))
})

## @cvxpy NONE
test_that("hermitian_wrap works with complex arg", {
  z <- Variable(c(3L, 3L), complex = TRUE)
  w <- hermitian_wrap(z)
  expect_true(is_hermitian(w))
})

# ── skew_symmetric_wrap ──────────────────────────────────────────────

## @cvxpy NONE
test_that("skew_symmetric_wrap validates real square", {
  x <- Variable(c(3L, 2L))
  expect_error(skew_symmetric_wrap(x), "square")

  z <- Variable(c(3L, 3L), complex = TRUE)
  expect_error(skew_symmetric_wrap(z), "real")
})

## @cvxpy NONE
test_that("skew_symmetric_wrap asserts skew_symmetric", {
  x <- Variable(c(3L, 3L))
  w <- skew_symmetric_wrap(x)
  expect_true(is_skew_symmetric(w))
})

# ── DCP propagation through complex atoms ────────────────────────────

## @cvxpy NONE
test_that("Re/Im/Conj propagate affine curvature", {
  x <- Variable(c(3L, 1L))
  expect_true(is_affine(Re(x)))
  expect_true(is_affine(Im(x)))
  expect_true(is_affine(Conj(x)))
})

## @cvxpy NONE
test_that("Re/Im/Conj propagate quadratic", {
  x <- Variable(c(3L, 1L))
  expect_true(is_quadratic(Re(x)))
  expect_true(is_quadratic(Im(x)))
  expect_true(is_quadratic(Conj(x)))
})

## @cvxpy NONE
test_that("Re/Im/Conj of constant are constant", {
  c1 <- Constant(1+2i)
  expect_true(is_constant(Re(c1)))
  expect_true(is_constant(Im(c1)))
  expect_true(is_constant(Conj(c1)))
})

# ── Complex propagation through atoms ───────────────────────────────

## @cvxpy NONE
test_that("AffAtom is_complex propagates from args", {
  x <- Variable(c(2L, 1L), complex = TRUE)
  y <- Variable(c(2L, 1L))
  ## Wrap is an AffAtom subclass
  wx <- Wrap(x)
  wy <- Wrap(y)
  expect_true(is_complex(wx))
  expect_false(is_complex(wy))
})

## @cvxpy NONE
test_that("AffAtom is_imag propagates from args", {
  x <- Variable(c(2L, 1L), imag = TRUE)
  y <- Variable(c(2L, 1L))
  wx <- Wrap(x)
  wy <- Wrap(y)
  expect_true(is_imag(wx))
  expect_false(is_imag(wy))
})
