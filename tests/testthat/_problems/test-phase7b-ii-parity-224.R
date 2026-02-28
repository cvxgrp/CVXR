# Extracted from test-phase7b-ii-parity.R:224

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(3)
gme <- geo_mean(x, approx = FALSE)
gma <- geo_mean(x, approx = TRUE)
expect_s3_class(gme, "CVXR::GeoMean")
expect_false(S7_inherits(gme, GeoMeanApprox))
