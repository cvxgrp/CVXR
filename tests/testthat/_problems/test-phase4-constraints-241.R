# Extracted from test-phase4-constraints.R:241

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
expect_equal(num_cones(constr), 1L)
