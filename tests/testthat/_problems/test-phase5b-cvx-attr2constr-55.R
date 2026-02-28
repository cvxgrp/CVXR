# Extracted from test-phase5b-cvx-attr2constr.R:55

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2, nonpos = TRUE)
p <- Problem(Maximize(sum_entries(x)), list(x >= -3))
ca <- CvxAttr2Constr()
