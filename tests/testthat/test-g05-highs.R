library(CVXR)
HIGHS_AVAILABLE  <- "HIGHS" %in% installed_solvers()

# Minimize:
#  x_0 +  x_1 + 3
# Subject to:
#               x_1 <=  7
#  5 <=  x_0 + 2x_1 <= 15
#  6 <= 3x_0 + 2x_1
#  0 <= x_0 <= 4
#  1 <= x_1
test_that("Test HIGHS LP example",{
  skip_on_cran()
  skip_if_not(HIGHS_AVAILABLE, "Skipping HIGHS test as it is not available.!")
  x <- Variable(2)
  A <- matrix(c(0,  1, 1,  2, -1, -2, -3, -2,
                1,  0, -1,  0, 0, -1),
              byrow = TRUE, ncol = 2)
  b <- c(7, 15, -5, -6, 4, 0, -1)
  
  obj <- matrix(c(1,1), nrow = 1) %*% x + 3
  p <- Problem(Minimize(obj), list(A %*% x <= b))
  res <- solve(p, solver = "HIGHS")
  res$getValue(x)
  expect_equal(res$status, "optimal")
  expect_equal(res$value, 5.75, tolerance = 1e-7)
  expect_equal(res$getValue(x), matrix(c(0.50, 2.25)), tolerance = 1e-7)
})

# Minimize:
#  -x_2 - 3x_3 + (1/2) * (2 x_1^2 - 2 x_1x_3 + 0.2 x_2^2 + 2 x_3^2)
# Subject to:
#  x_1 + x_3 <= 2
#  0 <= x
test_that("Test HIGHS QP example",{
  skip_on_cran()
  skip_if_not(HIGHS_AVAILABLE, "Skipping HIGHS test as it is not available.!")
  y <- Variable(3)
  L <- matrix(c(0, -1, -3), nrow = 1)
  Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
  A <- matrix(cbind(1, 0, 1), nrow = 1)
  obj <- L %*% y + 0.5 * quad_form(y, Q)
  constraints <- list(A %*% y <= 2, x >= 0)
  p <- Problem(Minimize(obj), constraints)
  res <- solve(p, solver = "HIGHS")
  res$getValue(y)
  expect_equal(res$status, "optimal")
  expect_equal(res$value, -5.25, tolerance = 1e-7)
  expect_equal(res$getValue(y), matrix(c(0.5, 5, 1.5)), tolerance = 1e-6)
})
