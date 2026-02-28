# Extracted from test-phase3e-affine-atoms.R:294

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(4)
expect_error(Convolve(x, y), "constant")
