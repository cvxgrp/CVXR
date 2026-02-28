# Extracted from test-phase4-constraints.R:209

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
expect_true(S7::S7_inherits(constr, PSD))
expect_true(S7::S7_inherits(constr, Cone))
