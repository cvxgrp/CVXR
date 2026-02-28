# Extracted from test-phase2c-operators.R:247

# test -------------------------------------------------------------------------
x <- Variable(1, nonneg = TRUE)
y <- Variable(1, nonpos = TRUE)
a <- AddExpression(list(x, y))
