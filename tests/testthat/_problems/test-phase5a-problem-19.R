# Extracted from test-phase5a-problem.R:19

# test -------------------------------------------------------------------------
x <- Variable(2)
obj <- Maximize(sum(x))
p <- Problem(obj, list(x <= 1))
expect_true(S7_inherits(p@objective, Maximize))
