x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test log problem", {
  # Log in objective
  obj <- Maximize(SumEntries(log(x)))
  constr <- list(x <= matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  
  # Log in constraint
  obj <- Minimize(SumEntries(x))
  constr <- list(log(x) >= 0, x <= matrix(c(1,1)))
  p <- Problem(obj, constr)
  
  # Index into log
  # obj <- Maximize(log(x)[2])
  # constr <- list(x <= matrix(c(1, exp(1))))
  # p <- Problem(obj, constr)
})

test_that("Test sigma max", {
  const <- Constant(cbind(c(1,2,3), c(4,5,6)))
  constr <- list(C == const)
  prob <- Problem(Minimize(Pnorm(C, 2)), constr)
})

test_that("Test a problem with exp", {
  for(n in c(5, 10, 25)) {
    x <- Variable(n)
    obj <- Minimize(SumEntries(exp(x)))
    p <- Problem(obj, list(SumEntries(x) == 1))
  }
})

test_that("Test a problem with log", {
  for(n in c(5, 10, 25)) {
    x <- Variable(n)
    obj <- Maximize(SumEntries(log(x)))
    p <- Problem(obj, list(SumEntries(x) == 1))
  }
})
