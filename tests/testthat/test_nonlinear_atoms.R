x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test log problems", {
  # Log in objective
  obj <- Maximize(SumEntries(log(x)))
  constr <- list(x <= as.matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  # result <- solve(p, solver = "CVXOPT")
  # expect_equal(result, 1, tolerance = 1e-6)
  # expect_equal(value(p, x), c(1, exp(1)), tolerance = 1e-6)
  
  # Log in constraint
  obj <- Minimize(SumEntries(x))
  constr <- list(log(x) >= 0, x <= as.matrix(c(1,1)))
  p <- Problem(obj, constr)
  # result <- solve(p, solver = "CVXOPT")
  # expect_equal(result, 2, tolerance = 1e-6)
  # expect_equal(value(p, x), c(1, 1), tolerance = 1e-6)

  # Index into log
  # obj <- Maximize(log(x)[1])
  # constr <- list(x <= as.matrix(c(1, exp(1))))
  # p <- Problem(obj, constr)
  # result <- solve(p, solver = "CVXOPT")
  # expect_equal(result, 1, tolerance = 1e-6)
  
  # Scalar log
  # obj <- Maximize(log(x[1]))
  # constr <- list(x <= as.matrix(c(1, exp(1))))
  # p <- Problem(obj, constr)
  # result <- solve(p, solver = "CVXOPT")
  # expect_equal(result, 1, tolerance = 1e-6)
})

test_that("Test a problem with entr", {
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Maximize(SumEntries(Entr(x)))
    p <- Problem(obj, list(SumEntries(x) == 1))
    # TODO: Finish this when all solvers are ready
  }
})

test_that("Test a problem with exp", {
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Minimize(SumEntries(exp(x)))
    p <- Problem(obj, list(SumEntries(x) == 1))
    # TODO: Finish this when all solvers are ready
  }
})

test_that("Test a problem with log", {
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Maximize(SumEntries(log(x)))
    p <- Problem(obj, list(SumEntries(x) == 1))
    # TODO: Finish this when all solvers are ready
  }
})

