# Extracted from test-phase3e-affine-atoms.R:285

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
a <- Constant(matrix(c(1, 2, 3), ncol = 1))
x <- Variable(4)
cv <- Convolve(a, x)
