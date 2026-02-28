# Extracted from test-phase5a-dcp2cone.R:33

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Minimize(sum(x)), list(x >= 0))
d2c <- Dcp2Cone()
