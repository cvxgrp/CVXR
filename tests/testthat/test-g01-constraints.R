context("test-g01-constraints")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

SOC <- CVXR:::SOC
save_value <- CVXR:::save_value

test_that("test the EqConstraint class", {
  skip_on_cran()
  constr <- x == z
  expect_equal(name(constr), "x == z")
  expect_equal(dim(constr), c(2,1))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))

  x <- save_value(x, 2)
  z <- save_value(z, 2)
  constr <- x == z
  expect_true(constr_value(constr))
  x <- save_value(x, 3)
  constr <- x == z
  expect_false(constr_value(constr))

  value(x) <- c(2,1)
  value(z) <- c(2,2)
  constr <- x == z
  expect_false(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,1)), tolerance = TOL)

  value(z) <- c(2,1)
  constr <- x == z
  expect_true(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,0)))
  expect_equal(residual(constr), matrix(c(0,0)))

  expect_error(x == y)
})

test_that("test the LeqConstraint class", {
  skip_on_cran()
  constr <- x <= z
  expect_equal(name(constr), "x <= z")
  expect_equal(dim(constr), c(2,1))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  x <- save_value(x, 1)
  z <- save_value(z, 2)
  constr <- x <= z
  expect_true(constr_value(constr))
  x <- save_value(x, 3)
  constr <- x <= z
  expect_false(constr_value(constr))

  value(x) <- c(2,1)
  value(z) <- c(2,0)
  constr <- x <= z
  expect_false(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,1)), tolerance = TOL)

  value(z) <- c(2,2)
  constr <- x <= z
  expect_true(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,0)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,0)), tolerance = TOL)

  expect_error(x <= y)
})

test_that("Test the PSD constraint %>>%", {
  skip_on_cran()
  constr <- A %>>% B
  expect_equal(name(constr), "A + -B >> 0")
  expect_equal(dim(constr), c(2,2))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  A <- save_value(A, rbind(c(2,-1), c(1,2)))
  B <- save_value(B, rbind(c(1,0), c(0,1)))
  constr <- A %>>% B
  expect_true(constr_value(constr))
  expect_equal(violation(constr), 0, tolerance = TOL)
  expect_equal(residual(constr), 0, tolerance = TOL)

  B <- save_value(B, rbind(c(3,0), c(0,3)))
  constr <- A %>>% B
  expect_false(constr_value(constr))
  expect_equal(violation(constr), 1, tolerance = TOL)
  expect_equal(residual(constr), 1, tolerance = TOL)

  expect_error(x %>>% 0, "Non-square matrix in positive definite constraint.")
})

test_that("Test the PSD constraint %<<%", {
  skip_on_cran()
  constr <- A %<<% B
  expect_equal(name(constr), "B + -A >> 0")
  expect_equal(dim(constr), c(2,2))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  B <- save_value(B, rbind(c(2,-1), c(1,2)))
  A <- save_value(A, rbind(c(1,0), c(0,1)))
  constr <- A %<<% B
  expect_true(constr_value(constr))
  A <- save_value(A, rbind(c(3,0), c(0,3)))
  constr <- A %<<% B
  expect_false(constr_value(constr))

  expect_error(x %<<% 0, "Non-square matrix in positive definite constraint.")
})

test_that("test the >= operator", {
  skip_on_cran()
  constr <- z >= x
  expect_equal(name(constr), "x <= z")
  expect_equal(dim(constr), c(2,1))
  expect_error(y >= x)
})

test_that("test the SOC class", {
  skip_on_cran()
  exp <- x + z
  scalar_exp <- a + b
  constr <- SOC(scalar_exp, exp)
  expect_equal(size(constr), 3)
})
