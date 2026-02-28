# Extracted from test-phase7c-expcone-atoms.R:84

# test -------------------------------------------------------------------------
x <- Variable(2, nonneg = TRUE)
y <- Variable(2, nonneg = TRUE)
e <- rel_entr(x, y)
s <- sign_from_args(e)
