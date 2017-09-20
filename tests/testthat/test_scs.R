TOL <- 1e-3

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test log problem", {
  # Log in objective
  obj <- Maximize(sum(log(x)))
  constr <- list(x <= matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, exp(1))), tolerance = TOL)
  
  # Log in constraint
  obj <- Minimize(sum(x))
  constr <- list(log(x) >= 0, x <= matrix(c(1,1)))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)
  
  # Index into log
  obj <- Maximize(log(x)[2])
  constr <- list(x <= c(1, exp(1)))
  p <- Problem(obj, constr)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 1, tolerance = TOL)
})

test_that("Test sigma max", {
  const <- Constant(rbind(c(1,2), c(3,4), c(5,6)))
  constr <- list(C == const)
  prob <- Problem(Minimize(norm(C, "F")), constr)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, value(sqrt(sum(const^2))), tolerance = TOL)
  expect_equal(result$getValue(C), value(const), tolerance = TOL)
})

test_that("Test sdp variable", {
  const <- Constant(rbind(c(1,2,3), c(4,5,6), c(7,8,9)))
  X <- Semidef(3)
  prob <- Problem(Minimize(0), list(X == const))
  result <- solve(prob, verbose = TRUE, solver = "SCS")
  expect_equal(tolower(result$status), "infeasible")
})

test_that("Test warm starting", {
  x <- Variable(10)
  obj <- Minimize(sum(exp(x)))
  prob <- Problem(obj, list(sum(x) == 1))
  result <- solve(prob, solver = "SCS")
  expect_equal(solve(prob, solver = "SCS")$value, result$value)
  # expect_false(solve(prob, solver = "SCS", warm_start = TRUE, verbose = TRUE) == result)
})
