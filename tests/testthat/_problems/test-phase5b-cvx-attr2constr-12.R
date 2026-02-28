# Extracted from test-phase5b-cvx-attr2constr.R:12

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
ca <- CvxAttr2Constr()
