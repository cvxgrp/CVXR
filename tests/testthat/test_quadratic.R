test_that("test elementwise power", {
  x <- Variable(3)
  y <- Variable(3)
  expect_false(is_constant(x))
  expect_true(is_affine(x))
  expect_true(is_quadratic(x))
  
  s <- Power(t(x) %*% y, 0)
  expect_true(is_constant(s))
  expect_true(is_affine(s))
  expect_true(is_quadratic(s))
  
  t <- Power(x-y, 1)
  expect_false(is_constant(t))
  expect_true(is_affine(t))
  expect_true(is_quadratic(t))
  
  u <- Power(x+2*y, 2)
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
  
  s <- t(x) %*% y
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

test_that("test QuadOverLin class", {
  x <- Variable(3, 5)
  y <- Variable(3, 5)
  z <- Variable()
  s <- QuadOverLin(x-y, z)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_false(is_quadratic(s))
  expect_true(is_dcp(s))
  
  t <- QuadOverLin(x+2*y, 5)
  expect_false(is_constant(t))
  expect_false(is_affine(t))
  expect_true(is_quadratic(t))
  expect_true(is_dcp(t))
})

test_that("test MatrixFrac class", {
  x <- Variable(5)
  M <- matrix(rnorm(25), nrow = 5, ncol = 5)
  P <- t(M) %*% M
  s <- MatrixFrac(x, P)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_true(is_dcp(s))
})

test_that("test quadratic form", {
  x <- Variable(5)
  P <- matrix(rnorm(25), nrow = 5, ncol = 5)
  q <- matrix(5, nrow = 5, ncol = 1)
  s <- t(x) %*% P %*% x + t(q) %*% x
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})

test_that("test SumSquares class", {
  X <- Variable(5, 4)
  P <- matrix(rnorm(15), nrow = 3, ncol = 5)
  Q <- matrix(rnorm(28), nrow = 4, ncol = 7)
  M <- matrix(rnorm(21), nrow = 3, ncol = 7)
  
  y <- P %*% X %*% Q + M
  expect_false(is_constant(y))
  expect_true(is_affine(y))
  expect_true(is_quadratic(y))
  expect_true(is_dcp(y))
  
  s <- SumSquares(y)
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_true(is_dcp(s))
  
  # Frobenius norm squared is indeed quadratic
  # but can't show quadraticity using recursive rules
  t <- norm(y, "fro")^2
  expect_false(is_constant(t))
  expect_false(is_affine(t))
  expect_false(is_quadratic(t))
  expect_true(is_dcp(t))
})

test_that("test indefinite quadratic", {
  x <- Variable()
  y <- Variable()
  z <- Variable()
  
  s <- y %*% z
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
  
  t <- (x+y)^2 - s - z %*% z
  expect_true(is_quadratic(t))
  expect_false(is_dcp(t))
})

test_that("test non-quadratic", {
  x <- Variable()
  y <- Variable()
  z <- Variable()
  
  expect_error(is_quadratic(x %*% y %*% z))
  
  s <- MaxEntries(VStack(x^2, Power(y, 2), z))
  expect_false(is_quadratic(s))
})

test_that("test affine product", {
  x <- Variable(3, 5)
  y <- Variable(5, 4)
  
  s <- x %*% y
  
  expect_false(is_constant(s))
  expect_false(is_affine(s))
  expect_true(is_quadratic(s))
  expect_false(is_dcp(s))
})
