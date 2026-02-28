# Extracted from test-phase3e-affine-atoms.R:165

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 1))
y <- Variable(c(3, 2))
h <- HStack(x, y)
