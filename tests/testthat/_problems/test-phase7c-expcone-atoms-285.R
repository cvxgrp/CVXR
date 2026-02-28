# Extracted from test-phase7c-expcone-atoms.R:285

# test -------------------------------------------------------------------------
x <- Variable(3, nonneg = TRUE)
e <- log_sum_exp(x)
s <- sign_from_args(e)
