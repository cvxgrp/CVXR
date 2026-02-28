# Extracted from test-phase1-expressions.R:838

# test -------------------------------------------------------------------------
x <- Variable(3)
expect_false(is_log_log_convex(x))
