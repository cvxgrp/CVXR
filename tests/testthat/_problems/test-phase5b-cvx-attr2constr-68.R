# Extracted from test-phase5b-cvx-attr2constr.R:68

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3, nonneg = TRUE)
y <- Variable(2)
expect_true("nonneg" %in% convex_attributes(list(x)))
