context("test-g01-dotsort")
TOL <- 1e-5

x <- Variable(5)

test_that("test sum k largest equivalence", {
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(1, 1, 1, 0)
  expr <- dotsort(x, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[(length(x_val)-3):length(x_val)]))
})

test_that("test sum k smallest equivalence", {
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(-1, -1, -1, 0)
  expr <- -dotsort(x, w)
  expect_true(is_concave(expr))
  expect_true(is_decr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[(length(x_val)-3):length(x_val)]))
})

test_that("test 1D", {
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(-1, 5, 2, 0, 5)
  expr <- dotsort(x, w)
  expect_true(is_convex(expr))
  expect_false(is_incr(expr, 1))
  expect_false(is_decr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(w)))
})

test_that("test 2D", {
  # X, W are matrix valued.
  xv <- Variable(5, 5)
  x_val <- matrix(0:24, nrow = 5, ncol = 5, byrow = TRUE)
  w <- matrix(0:3, nrow = 2, ncol = 2, byrow = TRUE)
  w_padded <- matrix(0, nrow = nrow(x_val), ncol = ncol(x_val))
  w_padded[1:nrow(w), 1:ncol(w)] <- w
  expr <- dotsort(x, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(xv == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(as.vector(x_val))*sort(as.vector(w_padded))))
})

test_that("test 0D", {
  # Scalar w.
  x_val <- c(1, 3, 2, -5, 0)
  w <- 1
  expr <- dotsort(x, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[length(x_val)]))
  
  # Scalar x, w.
  xv <- Variable()
  x_val <- c(1)
  w <- 1
  expr <- dotsort(xv, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(xv == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[length(x_val)]))
})




