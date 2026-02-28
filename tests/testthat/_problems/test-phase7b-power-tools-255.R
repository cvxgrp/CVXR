# Extracted from test-phase7b-power-tools.R:255

# test -------------------------------------------------------------------------
expect_error(fracify(c(1, -2, 3)), "nonnegative")
