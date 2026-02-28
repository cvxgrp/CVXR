# Extracted from test-phase2c-operators.R:310

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
a <- AddExpression(list(x, y))
copy <- expr_copy(a)
expect_equal(length(copy@args), 2L)
expect_true(S7_inherits(copy, AddExpression))
