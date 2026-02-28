# Extracted from test-phase5b-foundation.R:62

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
cons <- list(
    x >= 0,    # Inequality -> NonNeg from lower_ineq_to_nonneg
    x == 1     # Equality
  )
z <- Zero(x - Constant(1))
nn <- NonNeg(x)
result <- group_constraints(list(z, nn))
