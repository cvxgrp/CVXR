test_that("test 1D convolution", {
  n <- 3
  x <- Variable(n)
  f <- c(1, 2, 3)
  g <- c(0, 1, 0.5)
  f_conv_g <- c(0, 1, 2.5, 4, 1.5)
  expr <- Conv(f, g)
  expect_true(is_constant(expr))
  expect_equal(size(expr), c(5, 1))
  expect_equal(value(expr), f_conv_g)
  
  expr <- Conv(f, x)
  expect_true(is_affine(expr))
  expect_equal(size(expr), c(5, 1))
  
  # Matrix stuffing
  t <- Variable()
  prob <- Problem(Minimize(Norm(expr, 1)), list(x == g))
  result <- solve(prob)
  # expect_equal(result$optimal_value, sum(f_conv_g), tolerance = TOL)
  # expect_equal(value(expr), f_conv_g, tolerance = TOL)
})

test_that("test a problem with convolution", {
  N <- 5
  y <- matrix(rnorm(N), nrow = N, ncol = 1)
  h <- matrix(rnorm(2), nrow = 2, ncol = 1)
  x <- Variable(N)
  v <- Conv(h, x)
  obj <- Minimize(SumEntries(MulElemwise(y, v[1:N])))
  solve(Problem(obj, list()))
})
