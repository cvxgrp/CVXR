# Extracted from test-phase3c-axis-atoms.R:19

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
n <- Norm1(x, axis = 2L)
