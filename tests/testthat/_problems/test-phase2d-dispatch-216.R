# Extracted from test-phase2d-dispatch.R:216

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
z <- Variable(3)
w <- Variable(3)
result <- (x + y) + (z + w)
expect_true(S7_inherits(result, AddExpression))
