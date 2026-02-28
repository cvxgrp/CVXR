# Extracted from test-phase7b-power-tools.R:684

# test -------------------------------------------------------------------------
for (a in list(c(1, 2, 3), c(1, 1, 1, 1, 1), c(3, 7))) {
    result <- fracify(a)
    expect_true(check_dyad(result$w, result$w_dyad))
  }
