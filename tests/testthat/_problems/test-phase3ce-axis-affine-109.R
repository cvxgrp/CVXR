# Extracted from test-phase3ce-axis-affine.R:109

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
m <- MinEntries(x)
