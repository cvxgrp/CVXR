# Extracted from test-phase2g-objectives.R:157

# test -------------------------------------------------------------------------
x <- Variable(1, name = "x")
z <- Variable(1, name = "z")
exp <- x + z
obj <- Minimize(exp)
result <- canonical_form(obj)
