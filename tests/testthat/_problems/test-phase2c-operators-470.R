# Extracted from test-phase2c-operators.R:470

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3), PSD = TRUE)
y <- Variable(c(3, 3), PSD = TRUE)
m <- Multiply(x, y)
