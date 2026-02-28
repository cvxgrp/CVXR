# Extracted from test-cvxpy-atoms-parity.R:377

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
mat <- Constant(matrix(c(1, -1), nrow = 1))
expr <- SumEntries(mat %*% (Variable(2)^2))
