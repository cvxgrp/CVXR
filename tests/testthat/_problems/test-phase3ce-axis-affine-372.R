# Extracted from test-phase3ce-axis-affine.R:372

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
value(x) <- matrix(1:9, 3, 3)
u <- UpperTri(x)
