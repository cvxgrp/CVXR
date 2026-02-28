# Extracted from test-phase4-constraints.R:253

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
expect_equal(constr_size(constr), 9L)
