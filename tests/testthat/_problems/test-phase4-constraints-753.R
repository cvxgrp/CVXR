# Extracted from test-phase4-constraints.R:753

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3L, 1L))
constr <- NonPos(x)
expect_false(is_dgp(constr))
