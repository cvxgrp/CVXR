# Extracted from test-cvxpy-atoms-parity.R:854

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
for (axis in c(2L, 1L)) {
    x <- Variable(c(4L, 3L))
    x_val <- matrix(0:11, 4, 3)
    value(x) <- x_val

    expr <- Cumsum(x, axis = axis)
    ## R cumsum along axis
    if (axis == 2L) {
      target <- apply(x_val, 2, cumsum)
    } else {
      target <- t(apply(x_val, 1, cumsum))
    }
    expect_equal(as.matrix(value(expr)), target, tolerance = 1e-10,
                 info = paste("axis =", axis))
  }
