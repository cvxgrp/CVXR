# Extracted from test-phase5b-cvx-attr2constr.R:22

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2, nonneg = TRUE)
p <- Problem(Minimize(sum_entries(x)), list(x <= 3))
ca <- CvxAttr2Constr()
