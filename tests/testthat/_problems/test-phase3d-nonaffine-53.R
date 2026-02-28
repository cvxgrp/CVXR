# Extracted from test-phase3d-nonaffine.R:53

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(1, nonneg = TRUE)
value(x) <- c(3, 4, 0)
value(y) <- 5
q <- QuadOverLin(x, y)
