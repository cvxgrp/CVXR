# Extracted from test-cvxpy-atoms-parity.R:768

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
set.seed(42)
v <- rnorm(100)
x <- Constant(matrix(v, ncol = 1))
for (i in c(5, 10, 25, 50)) {
    expr <- SumLargest(x, k = i)
    prev_idx <- order(-v)[1:i]
    expect_equal(as.numeric(value(expr)), sum(v[prev_idx]), tolerance = 1e-6,
                 info = paste("k =", i))
  }
