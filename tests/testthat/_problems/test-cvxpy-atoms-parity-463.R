# Extracted from test-cvxpy-atoms-parity.R:463

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
entries <- list(x[1L, ], x[2L, ])
atom <- do.call(VStack, entries)
