# Extracted from test-phase5b-foundation.R:114

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3, nonneg = TRUE)
expect_true("nonneg" %in% convex_attributes(list(x)))
