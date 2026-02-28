## Tests for convenience atom wrappers (Category 1 missing atoms)

# ═══════════════════════════════════════════════════════════════════
# ptp (peak-to-peak)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("ptp: scalar result for vector", {
  x <- Variable(5)
  p <- ptp(x)
  expect_equal(p@shape, c(1L, 1L))
  expect_true(is_convex(p))
})

## @cvxpy NONE
test_that("ptp: axis=2 and axis=1", {
  x <- Variable(c(3, 4))
  p2 <- ptp(x, axis = 2L)
  p1 <- ptp(x, axis = 1L)
  expect_equal(p2@shape, c(1L, 4L))
  expect_equal(p1@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("ptp: solves correctly", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  ## Pin x = c(1, 2, 3) → ptp = max - min = 3 - 1 = 2
  prob <- Problem(Minimize(ptp(x)), list(x >= c(1, 2, 3), x <= c(1, 2, 3)))
  val <- psolve(prob)
  expect_equal(val, 2, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# cvxr_mean
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvxr_mean: full reduction is affine", {
  x <- Variable(c(2, 3))
  m <- cvxr_mean(x)
  expect_equal(m@shape, c(1L, 1L))
  expect_true(is_affine(m))
})

## @cvxpy NONE
test_that("cvxr_mean: axis=2", {
  x <- Variable(c(3, 2))
  m <- cvxr_mean(x, axis = 2L)
  expect_equal(m@shape, c(1L, 2L))
  expect_true(is_affine(m))
})

## @cvxpy NONE
test_that("cvxr_mean: axis=1", {
  x <- Variable(c(3, 2))
  m <- cvxr_mean(x, axis = 1L)
  expect_equal(m@shape, c(3L, 1L))
  expect_true(is_affine(m))
})

## @cvxpy NONE
test_that("cvxr_mean: solves correctly", {
  skip_if_not_installed("clarabel")
  x <- Variable(4)
  prob <- Problem(Minimize(cvxr_mean(x)), list(x >= c(1, 2, 3, 4)))
  val <- psolve(prob)
  expect_equal(val, 2.5, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# cvxr_std
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvxr_std: full reduction is convex", {
  x <- Variable(4)
  s <- cvxr_std(x)
  expect_equal(s@shape, c(1L, 1L))
  expect_true(is_convex(s))
})

## @cvxpy NONE
test_that("cvxr_std: axis=0 returns correct shape", {
  x <- Variable(c(3, 2))
  s <- cvxr_std(x, axis = 2L)
  expect_equal(s@shape, c(1L, 2L))
})

## @cvxpy NONE
test_that("cvxr_std: minimize std finds equal values", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob <- Problem(Minimize(cvxr_std(x)), list(sum(x) == 6, x >= 0))
  psolve(prob)
  xval <- as.numeric(value(x))
  ## All should be 2 (equal values minimize std)
  expect_equal(xval, c(2, 2, 2), tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# cvxr_var
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvxr_var: full reduction is convex", {
  x <- Variable(4)
  v <- cvxr_var(x)
  expect_equal(v@shape, c(1L, 1L))
  expect_true(is_convex(v))
})

## @cvxpy NONE
test_that("cvxr_var: axis not supported yet", {
  x <- Variable(c(3, 2))
  expect_error(cvxr_var(x, axis = 2L), "does not yet support")
})

## @cvxpy NONE
test_that("cvxr_var: minimize var finds equal values", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob <- Problem(Minimize(cvxr_var(x)), list(sum(x) == 9, x >= 0))
  psolve(prob)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(3, 3, 3), tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# vdot / scalar_product
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("vdot: shape is scalar", {
  x <- Variable(5)
  y <- Variable(5)
  d <- vdot(x, y)
  expect_equal(d@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("vdot: with constant", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  a <- c(1, 2, 3)
  prob <- Problem(Maximize(vdot(a, x)), list(x >= 0, x <= 1))
  val <- psolve(prob)
  ## Maximize 1*x1 + 2*x2 + 3*x3 with x in [0,1] → all x=1, val=6
  expect_equal(val, 6, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("vdot: matrix inputs are flattened", {
  x <- Variable(c(2, 3))
  y <- Variable(c(2, 3))
  d <- vdot(x, y)
  expect_equal(d@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("scalar_product: alias for vdot", {
  x <- Variable(3)
  y <- Constant(c(1, 2, 3))
  expect_true(S7::S7_inherits(scalar_product(x, y), Expression))
})

# ═══════════════════════════════════════════════════════════════════
# cvxr_outer
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvxr_outer: shape is (n, m)", {
  x <- Variable(3)
  y <- Variable(4)
  o <- cvxr_outer(x, y)
  expect_equal(o@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("cvxr_outer: rejects matrix input", {
  x <- Variable(c(2, 3))
  y <- Variable(4)
  expect_error(cvxr_outer(x, y), "must be a vector")
})

## @cvxpy NONE
test_that("cvxr_outer: solves correctly", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  y <- c(1, 2, 3)
  ## Minimize sum(outer(x, y)) = sum(x %*% t(y)) = sum(x) * sum(y) = 6*sum(x)
  prob <- Problem(Minimize(sum(cvxr_outer(x, y))), list(x >= 1))
  val <- psolve(prob)
  expect_equal(val, 12, tolerance = 1e-4)  # 2*6 = 12
})
