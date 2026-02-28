# Extracted from test-phase3ce-axis-affine.R:185

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
value(x) <- c(3, 4, 0)
n <- Pnorm(x, p = 2)
