# Extracted from test-cvxpy-parity.R:605

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
P <- Constant(matrix(c(1, 2, 3, 4), 2, 2))
expect_error(QuadForm(x, P), "[Ss]ymmetric|[Hh]ermitian")
