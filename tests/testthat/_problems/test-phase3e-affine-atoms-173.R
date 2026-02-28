# Extracted from test-phase3e-affine-atoms.R:173

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1))
y <- Variable(c(2, 1))
h <- HStack(x, y)
