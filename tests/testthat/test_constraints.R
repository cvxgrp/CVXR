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
})

test_that("test the LeqConstraint class", {
  constr <- x <= z
  expect_equal(size(constr), c(2,1))
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
