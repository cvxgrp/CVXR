# Extracted from test-phase3ce-axis-affine.R:421

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
c_atom <- cumsum_axis(x)
expect_true(S7_inherits(c_atom, Cumsum))
