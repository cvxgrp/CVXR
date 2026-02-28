# Extracted from test-phase5b-foundation.R:89

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
ineq <- (x <= Constant(c(3, 4)))
nn <- lower_ineq_to_nonneg(ineq)
