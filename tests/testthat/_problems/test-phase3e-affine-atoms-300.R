# Extracted from test-phase3e-affine-atoms.R:300

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
a <- Constant(matrix(c(1, 1), ncol = 1))
x <- Variable(3)
cv <- Convolve(a, x)
