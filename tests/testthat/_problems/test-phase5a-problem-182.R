# Extracted from test-phase5a-problem.R:182

# test -------------------------------------------------------------------------
x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
p <- Problem(Minimize(sum(x)), list(y >= 0))
inv <- InverseData(p)
