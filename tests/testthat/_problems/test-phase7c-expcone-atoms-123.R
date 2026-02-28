# Extracted from test-phase7c-expcone-atoms.R:123

# test -------------------------------------------------------------------------
x <- Variable(2, nonneg = TRUE)
e <- log1p_atom(x)
s <- sign_from_args(e)
