# Extracted from test-phase3ce-axis-affine.R:102

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
m <- max_entries(x, axis = 1L)
expect_true(S7_inherits(m, MaxEntries))
