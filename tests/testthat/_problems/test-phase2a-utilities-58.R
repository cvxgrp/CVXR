# Extracted from test-phase2a-utilities.R:58

# test -------------------------------------------------------------------------
v <- Variable(1, nonneg = TRUE)
c1 <- Constant(3)
result <- sum_signs(list(v, c1))
