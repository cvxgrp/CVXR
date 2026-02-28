# Extracted from test-phase3ce-axis-affine.R:405

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
c_atom <- Cumsum(x)
