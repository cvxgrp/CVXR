# Extracted from test-phase5a-problem.R:155

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Maximize(sum(x)), list(x >= 0, x <= 1))
fo <- FlipObjective()
