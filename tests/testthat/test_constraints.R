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
  
  x <- save_value(x, 2)
  z <- save_value(z, 2)
  constr <- x == z
  expect_true(value(constr))
  x <- save_value(x, 3)
  constr <- x == z
  expect_false(!is.na(value(constr)) && value(constr))
  
  value(x) <- c(2,1)
  value(z) <- c(2,2)
  constr <- x == z
  expect_false(!is.na(value(constr)) && value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(value(residual(constr)), matrix(c(0,1)), tolerance = TOL)
  
  value(z) <- c(2,1)
  constr <- x == z
  expect_true(value(constr))
  expect_equal(violation(constr), matrix(c(0,0)))
  expect_equal(value(residual(constr)), matrix(c(0,0)))
  
  expect_error(x == y)
})

test_that("test the LeqConstraint class", {
  constr <- x <= z
  expect_equal(size(constr), c(2,1))
  
  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_true(is.na(value(constr)))
  x <- save_value(x, 1)
  z <- save_value(z, 2)
  constr <- x <= z
  expect_true(value(constr))
  x <- save_value(x, 3)
  constr <- x <= z
  expect_false(!is.na(value(constr)) && value(constr))
  
  value(x) <- c(2,1)
  value(z) <- c(2,0)
  constr <- x <= z
  expect_false(!is.na(value(constr)) && value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(value(residual(constr)), matrix(c(0,1)), tolerance = TOL)
  
  value(z) <- c(2,2)
  constr <- x <= z
  expect_true(value(constr))
  expect_equal(violation(constr), matrix(c(0,0)), tolerance = TOL)
  expect_equal(value(residual(constr)), matrix(c(0,0)), tolerance = TOL)
  
  expect_error(x <= y)
})

test_that("Test the PSD constraint %>>%", {
  constr <- A %>>% B
  expect_equal(size(constr), c(2,2))
  
  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_true(is.na(value(constr)))
  A <- save_value(A, rbind(c(2,-1), c(1,2)))
  B <- save_value(B, rbind(c(1,0), c(0,1)))
  constr <- A %>>% B
  expect_true(value(constr))
  expect_equal(violation(constr), 0, tolerance = TOL)
  expect_equal(value(residual(constr)), 0, tolerance = TOL)
  
  B <- save_value(B, rbind(c(3,0), c(0,3)))
  constr <- A %>>% B
  expect_false(!is.na(value(constr)) && value(constr))
  expect_equal(violation(constr), 1, tolerance = TOL)
  expect_equal(value(residual(constr)), 1, tolerance = TOL)
  
  expect_error(x %>>% y)
})

test_that("Test the PSD constraint %<<%", {
  constr <- A %<<% B
  expect_equal(size(constr), c(2,2))
  
  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_true(is.na(value(constr)))
  B <- save_value(B, rbind(c(2,-1), c(1,2)))
  A <- save_value(A, rbind(c(1,0), c(0,1)))
  constr <- A %<<% B
  expect_true(value(constr))
  A <- save_value(A, rbind(c(3,0), c(0,3)))
  constr <- A %<<% B
  expect_false(!is.na(value(constr)) && value(constr))
  
  expect_error(x %<<% y)
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
