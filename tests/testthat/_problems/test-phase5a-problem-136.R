# Extracted from test-phase5a-problem.R:136

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Maximize(x), list(x <= 5))
fo <- FlipObjective()
