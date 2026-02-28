# Extracted from test-cvxpy-atoms-parity.R:243

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr <- MaxEntries(Variable(2), axis = 2L, keepdims = TRUE)
