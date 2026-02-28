# Extracted from test-v030-audit-fixes.R:267

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(4, 4))
m <- matrix(1:16, nrow = 4, ncol = 4)
dm0 <- DiagMat(x, k = 0L)
nv0 <- numeric_value(dm0, list(m))
