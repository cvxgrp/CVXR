# Extracted from test-phase5a-canonicalizers.R:224

# test -------------------------------------------------------------------------
x <- Variable(2)
y <- Variable(2)
t <- Variable(2)
result <- powcone_constrs(t, list(x, y), 0.5)
