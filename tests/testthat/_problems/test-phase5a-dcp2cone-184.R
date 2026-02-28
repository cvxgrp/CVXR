# Extracted from test-phase5a-dcp2cone.R:184

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(norm_inf(x)), list(x >= -1, x <= 1))
d2c <- Dcp2Cone()
