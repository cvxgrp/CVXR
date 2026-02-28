# Extracted from test-phase5a-problem.R:10

# test -------------------------------------------------------------------------
x <- Variable(2)
obj <- Minimize(sum(x))
p <- Problem(obj, list(x >= 0))
expect_true(S7_inherits(p, Problem))
