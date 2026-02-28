# Extracted from test-phase3e-affine-atoms.R:198

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
y <- Variable(c(1, 3))
v <- VStack(x, y)
