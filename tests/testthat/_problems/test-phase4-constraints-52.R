# Extracted from test-phase4-constraints.R:52

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
constr <- SOC(Variable(1L), x)
expect_true(S7::S7_inherits(constr, Cone))
