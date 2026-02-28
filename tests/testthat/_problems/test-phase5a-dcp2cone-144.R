# Extracted from test-phase5a-dcp2cone.R:144

# test -------------------------------------------------------------------------
x <- Variable(1)
y <- Variable(1)
p <- Problem(Minimize(y), list(exp(x) <= y, x >= -2))
d2c <- Dcp2Cone()
