# Extracted from test-phase3e-affine-atoms.R:128

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
d <- DiagMat(x)
val <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), 3, 3)
result <- numeric_value(d, list(val))
