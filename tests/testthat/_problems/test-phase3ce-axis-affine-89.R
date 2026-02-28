# Extracted from test-phase3ce-axis-affine.R:89

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
value(x) <- matrix(c(1, 5, 3, 2, 4, 6), 2, 3)
m <- MaxEntries(x)
