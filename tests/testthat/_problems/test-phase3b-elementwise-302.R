# Extracted from test-phase3b-elementwise.R:302

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
expect_error(Maximum(x), "at least 2")
