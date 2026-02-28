# Extracted from test-v030-audit-fixes.R:197

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
dv <- DiagVec(x, k = 1L)
expect_equal(dv@shape, c(4L, 4L))
expect_equal(dv@k, 1L)
expect_false(is_symmetric(dv))
expect_false(is_hermitian(dv))
