# Extracted from test-phase5a-dcp2cone.R:47

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(exp(x)), list(x >= -1))
d2c <- Dcp2Cone()
