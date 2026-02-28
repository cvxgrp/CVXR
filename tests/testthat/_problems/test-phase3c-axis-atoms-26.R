# Extracted from test-phase3c-axis-atoms.R:26

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
n <- Norm1(x, axis = 1L)
