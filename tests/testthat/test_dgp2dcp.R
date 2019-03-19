test_that("test unconstrained monomial", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod <- x*y
  dgp <- Problem(Minimize(prod))
  dgp2dcp <- Dgp2Dcp(dgp)
  
  dcp <- reduce(dgp2dcp)
  expect_equal(class(dcp@objective@expr), "AddExpression")
  expect_equal(length(dcp@objective@expr@args), 2)
  expect_equal(class(dcp@objective@expr@args[[1]]), "Variable")
  expect_equal(class(dcp@objective@expr@args[[2]]), "Variable")
  opt <- solve(dcp)
  
  # DCP is solved in log-space, so it is unbounded below
  # (since the OPT for dgp is 0 + epsilon).
  expect_equal(opt$value, -Inf)
  expect_equal(opt$status, "unbounded")
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, dcp@solution))
  expect_equal(dgp@value, 0.0)
  expect_equal(dgp@status, "unbounded")
  opt <- solve(dgp, gp = TRUE)
  expect_equal(opt$value, 0.0)
  expect_equal(opt$status, "unbounded")
  
  dgp <- Problem(Maximize(prod))
  dgp2dcp <- Dgp2Dcp(dgp)
  dcp <- reduce(dgp2dcp)
  opt <- solve(dcp)
  expect_equal(opt$value, Inf)
  expect_equal(opt$status, "unbounded")
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, dcp@solution))
  expect_equal(dgp@value, Inf)
  expect_equal(dgp@status, "unbounded")
  opt <- solve(dgp, gp = TRUE)
  expect_equal(opt$value, Inf)
  expect_equal(opt$status, "unbounded")
})

test_that("test basic equality constraint", {
  x <- Variable(pos = TRUE)
  dgp <- Problem(Minimize(x), list(x == 1.0))
  dgp2dcp <- Dgp2Dcp(dgp)
  
  dcp <- reduce(dgp2dcp)
  expect_equal(class(dcp@objective@expr), "Variable")
  opt <- solve(dcp)
  expect_equal(opt$value, 0.0, tolerance = TOL)
  expect_equal(value(variables(dcp)[[1]]), 0.0, tolerance = TOL)
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, dcp@solution))
  expect_equal(value(dgp), 1.0, tolerance = TOL)
  expect_equal(value(x), 1.0, tolerance = TOL)
  result <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
})
