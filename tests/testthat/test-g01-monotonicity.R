context("test-g01-monotonicity")

test_that("Test application of DCP composition rules to determine curvature", {
  expr <- 1 + exp(Variable())
  expect_equal(curvature(expr), CONVEX)
  
  expr <- Parameter() * NonNegative()
  expect_equal(curvature(expr), AFFINE)
  
  f <- function(x) { x^2 + x^0.5 }
  expr <- f(Constant(2))
  expect_equal(curvature(expr), CONSTANT)
  
  expr <- exp(Variable())^2
  expect_equal(curvature(expr), CONVEX)
  
  expr <- 1 - sqrt(Variable())
  expect_equal(curvature(expr), CONVEX)
  
  expr <- log(sqrt(Variable()))
  expect_equal(curvature(expr), CONCAVE)
  
  expr <- -(exp(Variable()))^2
  expect_equal(curvature(expr), CONCAVE)
  
  expr <- log(exp(Variable()))
  expect_false(is_dcp(expr))
  
  expr <- entr(NonNegative())
  expect_equal(curvature(expr), CONCAVE)
  
  expr <- ((Variable()^2)^0.5)^0
  expect_equal(curvature(expr), CONSTANT)
})

test_that("Test DCP composition rules with signed monotonicity", {
  # Convex argument
  expr <- abs(1 + exp(Variable()))
  expect_equal(curvature(expr), CONVEX)
  
  expr <- abs(-entr(Variable()))
  expect_equal(curvature(expr), UNKNOWN)
  
  expr <- abs(-log(Variable()))
  expect_equal(curvature(expr), UNKNOWN)
  
  # Concave argument
  expr <- abs(log(Variable()))
  expect_equal(curvature(expr), UNKNOWN)
  
  expr <- abs(square(Variable()))
  expect_equal(curvature(expr), CONVEX)
  
  expr <- abs(entr(Variable()))
  expect_equal(curvature(expr), UNKNOWN)
  
  # Affine argument
  expr <- abs(NonNegative())
  expect_equal(curvature(expr), CONVEX)
  
  expr <- abs(-NonNegative())
  expect_equal(curvature(expr), CONVEX)
  
  expr <- abs(Variable())
  expect_equal(curvature(expr), CONVEX)
})
