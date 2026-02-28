# Extracted from test-phase7c-expcone-atoms.R:381

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
e <- logistic(x)
result <- logistic_canon(e, e@args)
