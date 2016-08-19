test_that("Test application of DCP composition rules to determine curvature", {
  expr <- 1 + exp(Variable())
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- Parameter() * NonNegative()
  expect_equal(curvature(expr), Curvature.AFFINE)
  
  f <- function(x) { x^2 + x^0.5 }
  expr <- f(Constant(2))
  expect_equal(curvature(expr), Curvature.CONSTANT)
  
  expr <- exp(Variable())^2
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- 1 - sqrt(Variable())
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- log(sqrt(Variable()))
  expect_equal(curvature(expr), Curvature.CONCAVE)
  
  expr <- -(exp(Variable()))^2
  expect_equal(curvature(expr), Curvature.CONCAVE)
  
  expr <- log(exp(Variable()))
  expect_false(is_dcp(expr))
  
  expr <- Entr(NonNegative())
  expect_equal(curvature(expr), Curvature.CONCAVE)
  
  expr <- ((Variable()^2)^0.5)^0
  expect_equal(curvature(expr), Curvature.CONSTANT)
})

test_that("Test DCP composition rules with signed monotonicity", {
  # Convex argument
  expr <- abs(1 + exp(Variable()))
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- abs(-Entr(Variable()))
  expect_equal(curvature(expr), Curvature.UNKNOWN)
  
  expr <- abs(-log(Variable()))
  expect_equal(curvature(expr), Curvature.UNKNOWN)
  
  # Concave argument
  expr <- abs(log(Variable()))
  expect_equal(curvature(expr), Curvature.UNKNOWN)
  
  expr <- abs(Square(Variable()))
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- abs(Entr(Variable()))
  expect_equal(curvature(expr), Curvature.UNKNOWN)
  
  # Affine argument
  expr <- abs(NonNegative())
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- abs(-NonNegative())
  expect_equal(curvature(expr), Curvature.CONVEX)
  
  expr <- abs(Variable())
  expect_equal(curvature(expr), Curvature.CONVEX)
})
