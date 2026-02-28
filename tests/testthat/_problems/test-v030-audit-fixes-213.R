# Extracted from test-v030-audit-fixes.R:213

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
vals <- matrix(c(1, 2, 3), ncol = 1)
dv0 <- DiagVec(x, k = 0L)
nv0 <- numeric_value(dv0, list(vals))
