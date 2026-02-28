# Extracted from test-phase7b-power-tools.R:523

# test -------------------------------------------------------------------------
r <- pow_high(2)
expect_equal(sum(r$w), as.bigq(1L))
