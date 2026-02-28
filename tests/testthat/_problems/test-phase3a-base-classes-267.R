# Extracted from test-phase3a-base-classes.R:267

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
const <- create_const(matrix(1, 2, 2), c(2L, 2L))
