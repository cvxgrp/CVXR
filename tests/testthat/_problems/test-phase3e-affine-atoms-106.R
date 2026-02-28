# Extracted from test-phase3e-affine-atoms.R:106

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
d <- DiagVec(x)
val <- matrix(c(1, 2, 3), ncol = 1)
result <- numeric_value(d, list(val))
