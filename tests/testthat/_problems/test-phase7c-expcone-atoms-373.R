# Extracted from test-phase7c-expcone-atoms.R:373

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
e <- rel_entr(x, y)
result <- rel_entr_canon(e, e@args)
