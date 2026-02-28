# Extracted from test-phase2g-objectives.R:168

# test -------------------------------------------------------------------------
x <- Variable(1, name = "x")
z <- Variable(1, name = "z")
exp <- x + z
obj <- Maximize(exp)
result <- canonical_form(obj)
