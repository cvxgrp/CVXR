# Extracted from test-phase7c-expcone-atoms.R:23

# test -------------------------------------------------------------------------
x <- Variable(2, nonneg = TRUE)
y <- Variable(2, nonneg = TRUE)
e <- kl_div(x, y)
s <- sign_from_args(e)
