# Extracted from test-phase3e-affine-atoms.R:244

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2))
y <- Variable(c(2, 2))
expect_error(Kron(x, y), "constant")
