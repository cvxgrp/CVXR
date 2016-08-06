test_that("test 1D convolution", {
  n <- 3
  x <- Variable(n)
  f <- as.matrix(c(1, 2, 3))
  g <- as.matrix(c(0, 1, 0.5))
  f_conv_g <- as.matrix(c(0, 1, 2.5, 4, 1.5))
  expr <- Conv(f, g)
  expect_true(is_constant(expr))
  expect_equal(size(expr), c(5, 1))
  
  expr <- Conv(f, x)
  expect_true(is_affine(expr))
  expect_equal(size(expr), c(5, 1))
  
  t <- Variable()
  prob <- Problem(Minimize(Pnorm(expr, 1)), list(x == g))
})

test_that("test a problem with convolution", {
  N <- 5
  y <- matrix(rnorm(N), nrow = N, ncol = 1)
  h <- matrix(rnorm(2), nrow = 2, ncol = 1)
  x <- Variable(N)
  v <- Conv(h, x)
  # obj <- Minimize(SumEntries(MulElemwise(y, v[1:N])))
})
