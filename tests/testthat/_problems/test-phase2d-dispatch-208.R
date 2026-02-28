# Extracted from test-phase2d-dispatch.R:208

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
z <- Variable(3)
result <- x + (y + z)
expect_true(S7_inherits(result, AddExpression))
