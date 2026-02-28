# Extracted from test-phase5b-foundation.R:97

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
np <- NonPos(x)
nn <- nonpos2nonneg(np)
