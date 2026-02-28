# Extracted from test-phase5a-dcp2cone.R:107

# test -------------------------------------------------------------------------
x <- Variable(1)
y <- Variable(1)
p <- Problem(Minimize(exp(x) + abs(y)), list(x >= 0, y >= -1))
d2c <- Dcp2Cone()
