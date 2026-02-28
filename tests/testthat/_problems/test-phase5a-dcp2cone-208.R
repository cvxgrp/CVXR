# Extracted from test-phase5a-dcp2cone.R:208

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Minimize(sum(x)), list(x >= 0))
fo <- FlipObjective()
