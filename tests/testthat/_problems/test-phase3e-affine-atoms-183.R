# Extracted from test-phase3e-affine-atoms.R:183

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1))
y <- Variable(c(3, 1))
expect_error(HStack(x, y), "row")
