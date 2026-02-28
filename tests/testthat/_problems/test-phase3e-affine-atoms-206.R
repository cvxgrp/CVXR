# Extracted from test-phase3e-affine-atoms.R:206

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(1, 2))
y <- Variable(c(1, 2))
v <- VStack(x, y)
