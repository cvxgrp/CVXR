# Extracted from test-phase5b-foundation.R:105

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
z <- Zero(x - Constant(c(1, 2)))
expect_true(are_args_affine(list(z)))
