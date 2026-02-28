# Extracted from test-cvxpy-parity.R:611

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
P <- Constant(diag(2))
expect_error(QuadForm(x, P), "rows|match")
