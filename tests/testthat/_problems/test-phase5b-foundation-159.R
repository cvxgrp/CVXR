# Extracted from test-phase5b-foundation.R:159

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(SOLUTION_PRESENT, c("optimal", "optimal_inaccurate", "user_limit"))
