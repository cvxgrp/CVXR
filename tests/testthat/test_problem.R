a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the variables method", {
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars_ <- variables(p)
  ref <- list(a, x, b, A)
  mapply(function(v, r) { expect_equal(v, r) }, vars_, ref)
})

test_that("test the parameters method", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(p)
  ref <- c(p1, p2, p3)
  mapply(function(p, r) { expect_equal(p, r) }, params, ref)
})

test_that("test the constants method", {
  c1 <- matrix(randn(2), nrow = 1, ncol = 2)
  c2 <- matrix(randn(2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(c1*x), list(x >= c2))
  # constants_ <- constants(p)
  # ref <- list(as.character(c1), as.character(c2))
  # mapply(function(c, r) { expect_equal(c, r) }, constants_, ref)
})

test_that("test the is_dcp method", {
  p <- Problem(Minimize(NormInf(a)))
  expect_true(is_dcp(p))
  
  p <- Problem(Maximize(NormInf(a)))
  expect_false(is_dcp(p))
  # expect_error(solve(p))
  # solve(p, ignore_dcp = TRUE)
})

test_that("test problems involving variables with the same name", {
  var <- Variable(name = "a")
  p <- Problem(Maximize(a + var), list(var == 2 + a, var <= 3))
  # result <- solve(p)
  # expect_equal(result, 4.0, tolerance = TOL)
  # expect_equal(result$a, 1, tolerance = TOL)
  # expect_equal(result$var, 3, tolerance = TOL)
})

test_that("test adding problems", {
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob_minimize <- prob1 + prob2
  expect_equal(length(prob_minimize@constraints), 3)
  # result <- solve(prob_minimize)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)

  prob3 <- Problem(Maximize(a), list(b <= 1))
  prob4 <- Problem(Maximize(2*b), list(a <= 2))
  prob_maximize <- prob3 + prob4
  expect_equal(length(prob_maximize@constraints), 2)
  # result <- solve(prob_maximize)
  # expect_equal(result$optimal_value, 4, tolerance = TOL)
  
  # Test using the sum function
  prob5 <- Problem(Minimize(3*a))
  prob_sum <- Reduce("+", list(prob1, prob2, prob5))
  expect_equal(length(prob_sum@constraints), 3)
  # result <- solve(prob_sum)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  prob_sum <- Reduce("+", list(prob1))
  expect_equal(length(prob_sum@constraints), 1)
  
  # Test Minimize + Maximize
  expect_error(prob_bad_sum <- prob1 + prob3)
})

test_that("test problem multiplication by scalar", {
  prob1 <- Problem(Minimize(a^2), list(a >= 2))
  # answer <- solve(prob1)
  factors <- c(0, 1, 2.3, -4.321)
  for(f in factors) {
    # expect_equal(solve(f * prob1)$optimal_value, f * answer, tolerance = TOL)
    # expect_equal(solve(prob1 * f)$optimal_value, f * answer, tolerance = TOL)
  }
})

test_that("test problem linear combinations", {
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob3 <- Problem(Maximize(-(b + a)^2), list(b >= 3))
  
  # Simple addition and multiplication
  combo1 <- prob1 + 2 * prob2
  combo1_ref <- Problem(Minimize(a + 4*b), list(a >= b, a >= 1, b >= 2))
  # expect_equal(solve(combo1), solve(combo1_ref))
  
  # Division and subtraction
  combo2 <- prob1 - prob3/2
  combo2_ref <- Problem(Minimize(a + (b + a)^2), list(b >= 3, a >= b))
  # expect_equal(solve(combo2), solve(combo2_ref))
  
  # Multiplication with 0 (prob2's constraints should still hold)
  combo3 <- prob1 + 0*prob2 - 3*prob3
  combo3_ref <- Problem(Minimize(a + 3*(b + a)^2), list(a >= b, a >= 1, b >= 3))
  # expect_equal(solve(combo3), solve(combo3_ref))
})

test_that("test problems with parameters", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  p1@value <- 2
  p2@value <- -matrix(1, nrow = 3, ncol = 1)
  p3@value <- matrix(1, nrow = 4, ncol = 4)
  result <- solve(p)
  expect_equal(result$optimal_value, -6, tolerance = TOL)
  
  p1@value <- NA
  expect_error(solve(p))
})

test_that("test problems with NormInf", {
  # Constant argument
  p <- Problem(Minimize(NormInf(-2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  # Scalar arguments
  p <- Problem(Minimize(NormInf(a)), list(a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  
  p <- Problem(Minimize(3*NormInf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 3, tolerance = TOL)
  # expect_equal(result$a + 2*result$b, 0, tolerance = TOL)
  # expect_equal(result$c, 3, tolerance = TOL)

  # Maximize
  p <- Problem(Maximize(-NormInf(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, -2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Vector arguments
  p <- Problem(Minimize(NormInf(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  # expect_equal(result$x[2] - result$z[2], 7, tolerance = TOL)
})

test_that("Test problems with Norm1", {
  # Constant argument
  p <- Problem(Minimize(Norm1(-2)))
  # result <- solve(p)
  # expect_equal(result, 2, 1e-6)
  
  # Scalar arguments
  p <- Problem(Minimize(Norm1(a)), list(a <= -2))
  # result <- solve(p)
  
  # Maximize
  p <- Problem(Maximize(-Norm1(a)), list(a <= -2))
  # result <- solve(p)
  
  # Vector arguments
  p <- Problem(Minimize(Norm1(x - z) + 5), list(x >= matrix(c(2,3)), z <= matrix(c(-1,-4))))
  # result <- solve(p)
})

test_that("Test problems with Norm2", {
  # Constant argument
  p <- Problem(Minimize(Norm2(-2)))
  # result <- solve(p)
  # expect_equal(result, 2, 1e-6)
  
  # Scalar arguments
  p <- Problem(Minimize(Norm2(a)), list(a <= -2))
  # result <- solve(p)
  
  # Maximize
  p <- Problem(Maximize(-Norm2(a)), list(a <= -2))
  # result <- solve(p)
  
  # Vector arguments
  p <- Problem(Minimize(Norm2(x - z) + 5), list(x <= matrix(c(2,3)), z <= matrix(c(-1,-4))))
  # result <- solve(p)
  
  # Row arguments
  p <- Problem(Minimize(Norm2(t(x - z)) + 5), list(x >= matrix(c(2,3)), z <= matrix(c(-1,-4))))
  # result <- solve(p)
})

test_that("Test problems with abs", {
  p <- Problem(Minimize(SumEntries(abs(A))), list(-2 >= A))
  # result <- solve(p)
})
