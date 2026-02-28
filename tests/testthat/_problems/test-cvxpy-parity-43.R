# Extracted from test-cvxpy-parity.R:43

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
cvx <- Variable(1)^2
ccv <- Variable(1)^0.5
aff <- Variable(1)
const <- Constant(5)
expect_equal(expr_curvature(const - cvx), CONCAVE)
