# Extracted from test-cvxpy-atoms-parity.R:268

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr <- MinEntries(Variable(2), axis = 2L)
