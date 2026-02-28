# Extracted from test-phase5a-dcp2cone.R:16

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Maximize(sum(x)), list(x <= 1))
d2c <- Dcp2Cone()
