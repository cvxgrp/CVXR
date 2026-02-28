# Extracted from test-v030-audit-fixes.R:315

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
m <- matrix(1:9, nrow = 3, byrow = TRUE)
ut <- upper_tri(x)
nv <- numeric_value(ut, list(m))
