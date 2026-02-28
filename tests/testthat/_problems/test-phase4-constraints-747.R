# Extracted from test-phase4-constraints.R:747

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3L, 1L))
constr <- NonPos(x)
expect_equal(constr_size(constr), 3L)
