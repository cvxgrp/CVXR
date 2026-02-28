# Extracted from test-phase3e-affine-atoms.R:21

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
s <- SumEntries(x, axis = 2L)
