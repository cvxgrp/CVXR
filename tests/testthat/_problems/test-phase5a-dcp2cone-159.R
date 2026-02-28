# Extracted from test-phase5a-dcp2cone.R:159

# test -------------------------------------------------------------------------
x <- Variable(2)
cons <- list(x >= 0, x <= 5)
p <- Problem(Minimize(sum(x)), cons)
d2c <- Dcp2Cone()
