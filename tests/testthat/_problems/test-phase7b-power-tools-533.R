# Extracted from test-phase7b-power-tools.R:533

# test -------------------------------------------------------------------------
r <- pow_high(2, approx = FALSE)
expect_equal(sum(r$w), 1)
