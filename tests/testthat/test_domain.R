TOL <- 1e-6
a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")
z <- Variable(3, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test domain for partial minimization/maximization problems", {
  for(obj in list(Minimize(a^-1), Maximize(Log(a)))) {
    prob <- Problem(obj, list(x + a >= c(5,8)))
    # Optimize over nothing
    expr <- partial_optimize(prob, dont_opt_vars = list(x, a))
    dom <- domain(expr)
    constr <- list(a >= -100, x >= 0)
    prob <- Problem(Minimize(SumEntries(x + a)), c(dom, constr))
    result <- solve(prob)
    expect_equal(result$optimal_value, 13, tolerance = TOL)
    expect_true(result$a >= 0)
    expect_true(all(value(x + a - c(5, 8)) >= -1e-3))
    
    # Optimize over x
    expr <- partial_optimize(prob, opt_vars = list(x))
    dom <- domain(expr)
    constr <- list(a >= -100, x >= 0)
    prob <- Problem(Minimize(SumEntries(x + a)), c(dom, constr))
    result <- solve(prob)
    expect_equal(result$optimal_value, 0, tolerance = TOL)
    expect_true(result$a >= 0)
    expect_equal(result$x, c(0,0))
    
    # Optimize over x and a
    expr <- partial_optimize(prob, opt_vars = list(x, a))
    dom <- domain(expr)
    constr <- list(a >= -100, x >= 0)
    prob <- Problem(Minimize(SumEntries(x + a)), c(dom, constr))
    result <- solve(prob)
    expect_equal(result$a, -100, tolerance = TOL)
    expect_equal(result$x, c(0,0), tolerance = TOL)
  }
})

test_that("Test domain for GeoMean", {
  dom <- domain(GeoMean(x))
  prob <- Problem(Minimize(SumEntries(x)), dom)
  result <- solve(prob)
  expect_equal(result$optimal_value, 0, tolerance = TOL)
  
  # No special case for only one weight
  dom <- domain(GeoMean(x, c(0,2)))
  dom <- c(dom, x >= -1)
  prob <- Problem(Minimize(SumEntries(x)), dom)
  result <- solve(prob)
  expect_equal(result$x, c(-1,0), tolerance = TOL)
  
  dom <- domain(GeoMean(z, c(0,1,1)))
  dom <- c(dom, z >= -1)
  prob <- Problem(Minimize(SumEntries(z)), dom)
  result <- solve(prob)
  expect_equal(result$z, c(-1,0,0))
})

test_that("Test domain for QuadOverLin", {
  dom <- domain(QuadOverLin(x, a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$a, 0, tolerance = TOL)
})

test_that("Test domain for LambdaMax", {
  dom <- domain(LambdaMax(A))
  A0 <- rbind(c(1,2), c(3,4))
  result <- solve(Problem(Minimize(Norm2(A-A0)), dom))
  expect_equal(result$A, rbind(c(1,2.5), c(2.5,4)), tolerance = TOL)
})

test_that("Test domain for Pnorm", {
  dom <- domain(Pnorm(a, -0.5))
  prob <- Problem(Minimize(a), dom)
  result <- solve(prob)
  expect_equal(result$optimal_value, 0, tolerance = TOL)
})

test_that("Test domain for Log", {
  dom  <- domain(Log(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$a, 0, tolerance = TOL)
})

test_that("Test domain for Log1p", {
  dom <- domain(Log1p(a))
  result <- Problem(Minimize(a), dom)
  expect_equal(result$a, -1, tolerance = TOL)
})

test_that("Test domain for Entr", {
  dom <- domain(Entr(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$a, 0, tolerance = TOL)
})

test_that("Test domain for KLDiv", {
  b <- Variable()
  dom <- domain(KLDiv(a, b))
  result <- solve(Problem(Minimize(a + b), dom))
  expect_equal(result$a, 0, tolerance = TOL)
  expect_equal(result$b, 0, tolerance = TOL)
})

test_that("Test domain for Power", {
  dom <- domain(Sqrt(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$a, 0, tolerance = TOL)
  
  dom <- domain(Square(a))
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$a, -100, tolerance = TOL)
  
  dom <- domain(a^-1)
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$a, 0, tolerance = TOL)
  
  dom <- domain(a^3)
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$a, 0, tolerance = TOL)
})

test_that("Test domain for LogDet", {
  dom <- domain(LogDet(A + diag(rep(1,2))))
  prob <- Problem(Minimize(SumEntries(Diag(A))), dom)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$optimal_value, -2, tolerance = 1e-3)
})

test_that("Test domain for MatrixFrac", {
  dom <- domain(MatrixFrac(x, A + diag(rep(1,2))))
  prob <- Problem(Minimize(SumEntries(Diag(A))), dom)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$optimal_value, -2, tolerance = 1e-3)
})