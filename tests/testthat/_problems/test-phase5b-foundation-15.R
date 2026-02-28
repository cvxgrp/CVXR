# Extracted from test-phase5b-foundation.R:15

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
sol <- Solution(
    status = "optimal", opt_val = 42.0,
    primal_vars = list("1" = matrix(1:3)),
    dual_vars = list("2" = c(0.5, 0.5)),
    attr = list(solve_time = 0.01)
  )
