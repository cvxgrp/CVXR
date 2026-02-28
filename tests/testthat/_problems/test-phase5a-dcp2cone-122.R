# Extracted from test-phase5a-dcp2cone.R:122

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Maximize(log(x)), list(x >= 1, x <= 10))
fo <- FlipObjective()
