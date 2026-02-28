# Extracted from test-phase3e-affine-atoms.R:233

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
a <- Constant(diag(2))
x <- Variable(c(2, 2))
k <- Kron(a, x)
