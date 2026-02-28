# Extracted from test-phase3ce-axis-affine.R:389

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- c(1, 1, 1)
x <- Variable(5)
c_expr <- Convolve(a, x)
