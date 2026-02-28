# Extracted from test-phase5a-dcp2cone.R:172

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(norm1(x)), list(x >= -1, x <= 1))
d2c <- Dcp2Cone()
