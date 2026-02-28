# Extracted from test-phase7b-ii-parity.R:54

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
pe <- power(x, 2.5, approx = FALSE)
pa <- power(x, 2.5, approx = TRUE)
expect_s3_class(pe, "CVXR::Power")
expect_false(S7_inherits(pe, PowerApprox))
