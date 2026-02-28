# Extracted from test-phase7c-expcone-atoms.R:339

# test -------------------------------------------------------------------------
x <- Variable(c(3L, 4L))
e <- total_variation(x)
expect_true(S7_inherits(e, Expression))
