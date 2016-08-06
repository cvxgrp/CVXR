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

test_that("test the is_dcp method", {
  p <- Problem(Minimize(NormInf(a)))
  expect_true(is_dcp(p))
  
  p <- Problem(Maximize(NormInf(a)))
  expect_false(is_dcp(p))
})

test_that("test adding problems", {
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob_minimize <- prob1 + prob2
  expect_equal(length(prob_minimize@constraints), 3)
  # expect_equal(solve(prob_minimize), 6, 1e-6)

  prob3 <- Problem(Maximize(a), list(b <= 1))
  prob4 <- Problem(Maximize(2*b), list(a <= 2))
  prob_maximize <- prob3 + prob4
  expect_equal(length(prob_maximize@constraints), 2)
  # expect_equal(solve(prob_maximize), 4, 1e-6)
  
  # Test using the sum function
  prob5 <- Problem(Minimize(3*a))
  prob_sum <- Reduce("+", list(prob1, prob2, prob5))
  expect_equal(length(prob_sum@constraints), 3)
  # expect_equal(solve(prob_sum), 12, 1e-6)
  prob_sum <- Reduce("+", list(prob1))
  expect_equal(length(prob_sum@constraints), 1)
  
  # Test Minimize + Maximize
  expect_error(prob_bad_sum <- prob1 + prob3)
})

test_that("test problems with NormInf", {
  # Constant argument
  p <- Problem(Minimize(NormInf(-2)))
  # result <- solve(p)
  # expect_equal(result, 2, 1e-6)
  
  # Scalar arguments
  p <- Problem(Minimize(NormInf(a)), list(a >= 2))
  # result <- solve(p)
  # expect_equal(result, 2, 1e-6)
  
  p <- Problem(Minimize(3*NormInf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
  # result <- solve(p)
  # expect_equal(result, 3, 1e-6)

  # Maximize
  p <- Problem(Maximize(-NormInf(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result, -2, 1e-6)
  
  # Vector arguments
  p <- Problem(Minimize(NormInf(x - z) + 5), list(x >= matrix(c(2,3)), z <= matrix(c(-1,-4))))
  # result <- solve(p)
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
