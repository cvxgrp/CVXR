# Extracted from test-phase7c-expcone-atoms.R:207

# test -------------------------------------------------------------------------
x <- Variable(2, nonneg = TRUE)
e <- xexp(x)
expect_true(is_atom_log_log_convex(e))
