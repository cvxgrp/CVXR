# Extracted from test-phase5b-foundation.R:81

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
eq <- (x == Constant(c(1, 2)))
z <- lower_equality(eq)
