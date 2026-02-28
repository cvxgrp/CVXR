# Extracted from test-phase5a-dcp2cone.R:196

# test -------------------------------------------------------------------------
x <- Variable(5)
p <- Problem(Minimize(sum_largest(x, 2)), list(x >= -1, x <= 1))
d2c <- Dcp2Cone()
