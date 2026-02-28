# Extracted from test-phase4-constraints.R:234

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
nm <- expr_name(constr)
