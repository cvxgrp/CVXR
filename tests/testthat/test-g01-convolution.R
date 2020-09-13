context("test-g01-convolution")
TOL <- 1e-6

test_that("test 1D convolution", {
  skip_on_cran()
  n <- 3
  x <- Variable(n)
  f <- c(1, 2, 3)
  g <- c(0, 1, 0.5)
  f_conv_g <- c(0, 1, 2.5, 4, 1.5)
  expr <- conv(f, g)
  expect_true(is_constant(expr))
  expect_equal(dim(expr), c(5, 1))
  expect_equal(value(expr), f_conv_g)

  expr <- conv(f, x)
  expect_true(is_affine(expr))
  expect_equal(dim(expr), c(5, 1))

  # Matrix stuffing
  t <- Variable()
  prob <- Problem(Minimize(p_norm(expr, 1)), list(x == g))
  result <- solve(prob)
  expect_equal(result$value, sum(f_conv_g), tolerance = 1e-3)
  expect_equal(result$getValue(expr), f_conv_g, tolerance = TOL)
})

test_that("test a problem with convolution", {
  skip_on_cran()
  N <- 5
  y <- matrix(stats::rnorm(N), nrow = N, ncol = 1)
  h <- matrix(stats::rnorm(2), nrow = 2, ncol = 1)
  x <- Variable(N)
  v <- conv(h, x)
  obj <- Minimize(sum(y * v[1:N]))
  result <- solve(Problem(obj, list()))
  expect_equal(result$status, "unbounded")
})
