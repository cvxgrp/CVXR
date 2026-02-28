# Extracted from test-phase7c-expcone-atoms.R:393

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
e <- log1p_atom(x)
result <- log1p_canon(e, e@args)
