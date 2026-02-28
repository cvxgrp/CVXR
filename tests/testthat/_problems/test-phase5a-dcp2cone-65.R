# Extracted from test-phase5a-dcp2cone.R:65

# test -------------------------------------------------------------------------
x <- Variable(1, nonneg = TRUE)
p <- Problem(Minimize(-log(x)), list(x >= 1))
d2c <- Dcp2Cone()
