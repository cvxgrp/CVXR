# Extracted from test-phase5a-dcp2cone.R:94

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(max_entries(x)), list(x >= -1, x <= 1))
d2c <- Dcp2Cone()
