# Extracted from test-phase7b-ii-parity.R:109

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
pne <- p_norm(x, p = 3, approx = FALSE)
pna <- p_norm(x, p = 3, approx = TRUE)
expect_s3_class(pne, "CVXR::Pnorm")
expect_false(S7_inherits(pne, PnormApprox))
