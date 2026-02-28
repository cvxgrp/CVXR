# Extracted from test-cvxpy-atoms-parity.R:1110

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
for (n in 3:7) {
    A <- upper_tri_to_full(n)
    ell <- (n * (n + 1L)) %/% 2L
    v <- seq_len(ell) - 1L  # 0-based like np.arange
    ## A @ v reshaped to (n, n) F-order should be symmetric
    Mv <- as.numeric(A %*% v)
    M <- matrix(Mv, n, n)  # F-order (column-major)
    expect_equal(M, t(M), tolerance = 1e-10, info = paste("n =", n))
  }
