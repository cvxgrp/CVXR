# Extracted from test-phase3e-affine-atoms.R:41

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
s <- SumEntries(x, axis = 2L)
