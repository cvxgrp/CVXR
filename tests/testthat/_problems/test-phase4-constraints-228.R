# Extracted from test-phase4-constraints.R:228

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
expect_false(is_dgp(constr))
