# Extracted from test-phase3b-elementwise.R:293

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
m <- Maximum(x, y)
