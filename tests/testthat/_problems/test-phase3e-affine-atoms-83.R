# Extracted from test-phase3e-affine-atoms.R:83

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
expect_error(Reshape(x, c(4, 2)), "does not match")
