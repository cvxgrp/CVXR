# Extracted from test-cvxpy-parity.R:66

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
cvx <- Variable(1)^2
ccv <- Variable(1)^0.5
aff <- Variable(1)
const <- Constant(5)
zero <- Constant(0)
neg <- Constant(-1)
pos <- Constant(1)
e_zero_cvx <- zero * cvx
expect_true(is_affine(e_zero_cvx))
expect_equal(expr_curvature(neg * cvx), CONCAVE)
