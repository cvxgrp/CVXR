# Extracted from test-phase3e-affine-atoms.R:225

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
a <- Constant(diag(2))
x <- Variable(c(3, 3))
k <- Kron(a, x)
