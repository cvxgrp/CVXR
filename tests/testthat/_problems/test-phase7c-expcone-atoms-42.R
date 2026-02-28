# Extracted from test-phase7c-expcone-atoms.R:42

# test -------------------------------------------------------------------------
x <- Variable(2)
y <- Variable(2)
e <- kl_div(x, y)
dom <- atom_domain(e)
