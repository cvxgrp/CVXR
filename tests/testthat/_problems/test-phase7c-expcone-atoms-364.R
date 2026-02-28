# Extracted from test-phase7c-expcone-atoms.R:364

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
e <- kl_div(x, y)
result <- kl_div_canon(e, e@args)
