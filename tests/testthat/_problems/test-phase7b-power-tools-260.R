# Extracted from test-phase7b-power-tools.R:260

# test -------------------------------------------------------------------------
for (a in list(c(1, 2, 3), c(1, 1), c(3, 5, 7, 11))) {
    result <- fracify(a)
    expect_true(is_weight(result$w))
    expect_true(is_dyad_weight(result$w_dyad))
    expect_true(check_dyad(result$w, result$w_dyad))
  }
