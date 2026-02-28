# Extracted from test-cvxpy-atoms-parity.R:652

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
X_psd <- Variable(c(2L, 2L), PSD = TRUE)
X_nsd <- Variable(c(2L, 2L), NSD = TRUE)
expect_true(is_nonneg(Trace(X_psd)))
