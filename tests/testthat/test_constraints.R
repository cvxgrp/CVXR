TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the EqConstraint class", {
  constr <- x == z
  expect_equal(size(constr), c(2,1))
  
  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_true(is.na(value(constr)))
  
  x@primal_value <- 2
  z@primal_value <- 2
  expect_true(value(constr))
  x@primal_value <- 3
  expect_false(value(constr))
  
  x@primal_value <- c(2,1)
  z@primal_value <- c(2,2)
  expect_false(value(constr))
  expect_equal(violation(constr), c(0,1), tolerance = TOL)
  expect_equal(value(constr@residual), c(0,1), tolerance = TOL)
  
  z@primal_value <- c(2,1)
  expect_true(value(constr))
  expect_equal(violation(constr), c(0,0))
  expect_equal(value(residual(constr)), c(0,0))
  
  expect_error(x == y)
})

test_that("test the LeqConstraint class", {
  constr <- x <= z
  expect_equal(size(constr), c(2,1))
  
  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_true(is.na(value(constr)))
  x@primal_value <- 1
  z@primal_value <- 2
  expect_true(value(constr))
  x@primal_value <- 3
  expect_false(value(constr))
  
  x@primal_value <- c(2,1)
  z@primal_value <- c(2,0)
  expect_false(value(constr))
  expect_equal(violation(constr), c(0,1), tolerance = TOL)
  expect_equal(value(residual(constr)), c(0,1), tolerance = TOL)
  
  z@primal_value <- c(2,2)
  expect_true(value(constr))
  expect_equal(violation(constr), c(0,0), tolerance = TOL)
  expect_equal(value(residual(constr)), c(0,0), tolerance = TOL)
  
  expect_error(x <= y)
})

test_that("Test the PSD constraint >>", {
  # constr <- A >> B
  # expect_equal(size(constr), c(2,2))
  
  # Test value and dual_value
  # expect_true(is.na(dual_value(constr)))
  # expect_true(is.na(value(constr)))
  # A@primal_value <- rbind(c(2,-1), c(1,2))
  # B@primal_value <- rbind(c(1,0), c(0,1))
  # expect_true(value(constr))
  # expect_equal(violation(constr), 0, tolerance = TOL)
  # expect_equal(value(residual(constr)), 0, tolerance = TOL)
  
  # B@primal_value <- rbind(c(3,0), c(0,3))
  # expect_false(value(constr))
  # expect_equal(violation(constr), 1, tolerance = TOL)
  # expect_equal(value(residual(constr)), 1, tolerance = TOL)
  
  # expect_error(x >> y)
})

test_that("Test the PSD constraint <<", {
  # constr <- A << B
  # expect_equal(size(constr), c(2,2))
  
  # Test value and dual_value
  # expect_true(is.na(dual_value(constr)))
  # expect_true(is.na(value(constr)))
  # B@primal_value <- rbind(c(2,-1), c(1,2))
  # A@primal_value <- rbind(c(1,0), c(0,1))
  # expect_true(value(constr))
  # A@primal_value <- rbind(c(3,0), c(0,3))
  # expect_false(value(constr))
  
  # expect_error(x << y)
})

test_that("test the < operator", {
  constr <- x < z
  expect_equal(size(constr), c(2,1))
  expect_error(x < y)
})

test_that("test the >= operator", {
  constr <- z >= x
  expect_equal(size(constr), c(2,1))
  expect_error(y >= x)
})

test_that("test the > operator", {
  constr <- z > x
  expect_equal(size(constr), c(2,1))
  expect_error(y > x)
})

test_that("test the SOC class", {
  exp <- x + z
  scalar_exp <- a + b
  constr <- SOC(scalar_exp, list(exp))
  expect_equal(size(constr), c(3,1))
})
