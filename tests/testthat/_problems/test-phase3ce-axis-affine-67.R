# Extracted from test-phase3ce-axis-affine.R:67

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
s <- sum_entries(x, axis = 1L)
expect_true(S7_inherits(s, SumEntries))
