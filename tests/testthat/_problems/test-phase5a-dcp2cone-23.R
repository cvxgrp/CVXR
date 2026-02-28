# Extracted from test-phase5a-dcp2cone.R:23

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(-exp(x)))
d2c <- Dcp2Cone()
