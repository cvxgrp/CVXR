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
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[(length(x_val)-3):length(x_val)]), tolerance = TOL)
})

test_that("test sum k smallest equivalence", {
  x_val <- c(1, 3, 2, -5, 0)
  w <- c(-1, -1, -1, 0)
  expr <- -dotsort(x, w)
  expect_true(is_concave(expr))
  expect_true(is_decr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[(length(x_val)-3):length(x_val)]), tolerance = TOL)
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
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(w)), tolerance = TOL)
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
  expect_equal(result$getValue(prob@objective), sum(sort(as.vector(x_val))*sort(as.vector(w_padded))), tolerance = TOL)
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
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[length(x_val)]), tolerance = TOL)
  
  # Scalar x, w.
  xv <- Variable()
  x_val <- c(1)
  w <- 1
  expr <- dotsort(xv, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(xv == x_val))
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)[length(x_val)]), tolerance = TOL)
})

test_that("test constant", {
  xc <- 0:24
  x_val <- matrix(0:24, nrow = 5, ncol = 5, byrow = TRUE)
  w <- matrix(0:3, nrow = 2, ncol = 2, byrow = TRUE)
  w_padded <- matrix(0, nrow = nrow(x_val), ncol = ncol(x_val))
  w_padded[1:nrow(w), 1:ncol(w)] <- w
  expr <- dotsort(x, w)
  expect_true(is_convex(expr))
  expect_true(is_incr(expr, 1))
  
  prob <- Problem(Minimize(expr), list())
  result <- solve(prob)
  expect_equal(result$getValue(prob@objective), sum(sort(as.vector(x_val))*sort(as.vector(w_padded))), tolerance = TOL)
})

test_that("test parameter", {
  x_val <- c(1, 3, 2, -5, 0)
  
  expect_true(is_incr(dotsort(x, Parameter(2, pos = TRUE)), 1))
  expect_true(is_incr(dotsort(x, Parameter(2, nonneg = TRUE)), 1))
  expect_false(is_incr(dotsort(x, Parameter(2, neg = TRUE)), 1))
  expect_true(is_decr(dotsort(x, Parameter(2, neg = TRUE)), 1))
  
  w_p <- Parameter(2, value = c(1, 0))
  expr <- dotsort(x, w_p)
  expect_false(is_incr(expr, 1))
  expect_false(is_decr(expr, 1))
  
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob, enforce_dpp = TRUE)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(c(1, 0, 0, 0, 0))), tolerance = TOL)
  value(w_p) <- c(-1, -1)
  # expr <- dotsort(x, w_p)
  # prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob, enforce_dpp = TRUE)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(c(-1, -1, 0, 0, 0))), tolerance = TOL)
  
  # Test parameter affine.
  w_p <- Parameter(2, value = c(1, 0))
  parameter_affine_expression <- 2*w_p
  expr <- dotsort(x, parameter_affine_expression)
  prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob, enforce_dpp = TRUE)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(c(2, 0, 0, 0, 0))), tolerance = TOL)
  value(w_p) <- c(-1, -1)
  # expr <- dotsort(x, w_p)
  # prob <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(prob, enforce_dpp = TRUE)
  expect_equal(result$getValue(prob@objective), sum(sort(x_val)*sort(c(-2, -2, 0, 0, 0))), tolerance = TOL)
  
  # Test constant + non-parameter affine.
  x_const <- c(1, 2, 3)
  p <- Parameter(value = 2)
  p_squared <- p^2
  expr <- dotsort(x_const, p_squared)
  problem <- Problem(Minimize(expr))
  result <- solve(problem, enforce_dpp = TRUE)
  expect_equal(result$getValue(expr), 2^2*3, tolerance = TOL)
  value(p) <- -1
  # p_squared <- p^2
  # expr <- dotsort(x_const, p_squared)
  # problem <- Problem(Minimize(expr))
  result <- solve(problem, enforce_dpp = TRUE)
  expect_equal(result$getValue(expr), (-1)^2*3)
  
  # Test non-parameter affine with implicit casting to constant.
  x_val <- c(1, 2, 3, 4, 5)
  p <- Parameter(value = 2)
  p_squared <- p^2
  expr <- dotsort(x, p_squared)
  problem <- Problem(Minimize(expr), list(x == x_val))
  result <- solve(problem)
  expect_equal(result$getValue(expr), 2^2*5, tolerance = TOL)
  value(p) <- -1
  # p_squared <- p^2
  # expr <- dotsort(x, p_squared)
  # problem <- Problem(Minimize(expr))
  result <- solve(problem)
  expect_equal(result$getValue(expr), (-1)^2*5, tolerance = TOL)
})

test_that("test list", {
  r <- matrix(c(2, 1, 0, -1, -1))
  w <- matrix(c(1.2, 1.1))
  expr <- dotsort(x, w)
  prob <- Problem(Maximize(t(r) %*% x), list(0 <= x, expr <= 1, sum(x) == 1))
  result <- solve(prob)
  expect_equal(result$getValue(expr), 1, tolerance = TOL)
  expect_equal(sum(result$getValue(x)[1:2] * w), 1, tolerance = TOL)
})

test_that("test composition", {
  r <- matrix(c(2, 1, 0, -1, -1))
  w <- matrix(c(0.7, 0.8))
  expr <- dotsort(x, w)
  prob <- Problem(Maximize(t(r) %*% x), list(0 <= x, expr <= 2, sum(x) == 1))
  result <- solve(prob)
  expect_equal(result$getValue(expr), 2, tolerance = TOL)
  expect_equal(sum(sort(exp(result$getValue(x))[(nrow(x) - 2):nrow(x)]) * sort(w)), 2, tolerance = TOL)
})

test_that("test non-fixed x", {
  r <- matrix(c(2, 1, 0, -1, -1))
  w <- matrix(c(1.2, 1.1))
  expr <- dotsort(x, w)
  prob <- Problem(Maximize(t(r) %*% x), list(0 <= x, expr <= 1, sum(x) == 1))
  result <- solve(prob)
  expect_equal(result$getValue(expr), 1, tolerance = TOL)
  expect_equal(sum(result$getValue(x)[1:2] * w), 1, tolerance = TOL)
  
  # Test non-fixed x - unordered w.
  r <- matrix(c(2, 1, 0, -1, -1))
  w <- matrix(c(1.2, 1.1, 1.3))
  expr <- dotsort(x, w)
  prob <- Problem(Maximize(t(r) %*% x), list(0 <= x, expr <= 1, sum(x) == 1))
  result <- solve(prob)
  expect_equal(result$getValue(expr), 1, tolerance = TOL)
  expect_equal(sum(result$getValue(x)[(nrow(x) - 3):nrow(x)] * w), 1, tolerance = TOL)
})

test_that("test exceptions", {
  # length(w) > length(x).
  expect_error(dotsort(x, c(1, 2, 3, 4, 5, 8)), 
               "The size of W must be less than or equal to the size of X", fixed = TRUE)
  
  # Two variable expressions.
  expect_error(dotsort(x, Variable(3)), 
               "The W argument must be constant", fixed = TRUE)
  
  # Swapped arguments.
  expect_error(dotsort(c(1, 2, 3), x), "
               The W argument must be constant", fixed = TRUE)
  
  # Non-DCP composition.
  expect_error(solve(Problem(Minimize(dotsort(abs(x), c(-1, 1))))), 
               "Problem does not follow DCP rules", fixed = TRUE)
  
  # Non-DPP composition.
  p <- Parameter(value = 2)
  p_squared <- p^2
  expect_error(solve(Problem(Minimize(dotsort(x, p_squared))), enforce_dpp = TRUE),
               "You are solving a parametrized problem that is not DPP", fixed = TRUE)
})
