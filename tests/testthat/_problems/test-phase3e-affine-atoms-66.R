# Extracted from test-phase3e-affine-atoms.R:66

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
r <- Reshape(x, c(-1, 6))
