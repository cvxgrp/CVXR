# Extracted from test-phase3ce-axis-affine.R:414

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(4)
value(x) <- c(1, 2, 3, 4)
c_atom <- Cumsum(x)
