# Extracted from test-phase5a-problem.R:146

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(x), list(x >= 0))
fo <- FlipObjective()
