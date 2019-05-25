context("test-g01-quadratic")

test_that("test elementwise power", {
  x <- Variable(3)
  y <- Variable(3)
  expect_false(is_constant(x))
  expect_true(is_affine(x))
  expect_true(is_quadratic(x))

  expect_warning(s <- power(t(x) %*% y, 0))
  expect_true(is_constant(s))
  expect_true(is_affine(s))
  expect_true(is_quadratic(s))

  t <- power(x-y, 1)
  expect_false(is_constant(t))
  expect_true(is_affine(t))
  expect_true(is_quadratic(t))

  u <- power(x+2*y, 2)
  expect_false(is_constant(u))
  expect_false(is_affine(u))
  expect_true(is_quadratic(u))
  expect_true(is_dcp(u))

  w <- (x+2*y)^2
  expect_false(is_constant(w))
  expect_false(is_affine(w))
  expect_true(is_quadratic(w))
  expect_true(is_dcp(w))
})

test_that("test matrix multiplication", {
  x <- Variable(3, 5)
  y <- Variable(3, 5)
  expect_false(is_constant(x))
  expect_true(is_affine(x))
  expect_true(is_quadratic(x))

  expect_warning(s <- t(x) %*% y)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

test_that("test quad_over_lin function", {
  x <- Variable(3, 5)
  y <- Variable(3, 5)
  z <- Variable()
  s <- quad_over_lin(x-y, z)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_false(is_quadratic(s))
  expect_true(is_dcp(s))

  t <- quad_over_lin(x+2*y, 5)
  expect_false(is_constant(t))
  expect_false(is_affine(t))
  expect_true(is_quadratic(t))
  expect_true(is_dcp(t))
})

test_that("test matrix_frac function", {
  x <- Variable(5)
  M <- diag(5)
  P <- t(M) %*% M
  s <- matrix_frac(x, P)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_true(is_dcp(s))
})

test_that("test quadratic form", {
  x <- Variable(5)
  P <- diag(5) - 2*matrix(1, nrow = 5, ncol = 5)
  q <- matrix(1, nrow = 5, ncol = 1)

  expect_warning(s <- t(x) %*% P %*% x + t(q) %*% x)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

test_that("test sum_squares function", {
  X <- Variable(5, 4)
  P <- matrix(1, nrow = 3, ncol = 5)
  Q <- matrix(1, nrow = 4, ncol = 7)
  M <- matrix(1, nrow = 3, ncol = 7)

  y <- P %*% X %*% Q + M
  expect_false(is_constant(y))
  expect_true(is_affine(y))
  expect_true(is_quadratic(y))
  expect_true(is_dcp(y))

  s <- sum_squares(y)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_true(is_dcp(s))

  # Frobenius norm squared is indeed quadratic
  # but can't show quadraticity using recursive rules
  t <- norm(y, "F")^2
  expect_false(is_constant(t))
  expect_false(is_affine(t))
  expect_false(is_quadratic(t))
  expect_true(is_dcp(t))
})

test_that("test indefinite quadratic", {
  x <- Variable()
  y <- Variable()
  z <- Variable()

  s <- y*z
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))

  t <- (x+y)^2 - s - z*z
  expect_true(is_quadratic(t))
  expect_false(is_dcp(t))
})

test_that("test non-quadratic", {
  x <- Variable()
  y <- Variable()
  z <- Variable()

  s <- max_entries(vstack(x, y, z))^2
  expect_false(is_quadratic(s))

  s <- max_entries(vstack(x^2, power(y, 2), z))
  expect_false(is_quadratic(s))
})

test_that("test affine product", {
  x <- Variable(3, 5)
  y <- Variable(5, 4)

  expect_warning(s <- x %*% y)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})
