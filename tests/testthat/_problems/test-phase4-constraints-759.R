# Extracted from test-phase4-constraints.R:759

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
constr <- NonPos(x)
save_dual_value(constr, matrix(c(1, 2), 2, 1))
