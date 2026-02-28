# Extracted from test-phase1-expressions.R:286

# test -------------------------------------------------------------------------
m <- matrix(1:6, 2, 3)
c1 <- Constant(m)
canon <- canonicalize(c1)
