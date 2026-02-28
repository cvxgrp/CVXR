# Extracted from test-phase2d-dispatch.R:68

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- +x
expect_true(S7_inherits(result, Expression))
