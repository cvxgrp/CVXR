# Extracted from test-phase3c-axis-atoms.R:123

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
m <- MaxEntries(x, axis = 2L)
