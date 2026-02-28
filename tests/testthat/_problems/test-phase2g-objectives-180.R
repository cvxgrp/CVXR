# Extracted from test-phase2g-objectives.R:180

# test -------------------------------------------------------------------------
x <- Variable(1, name = "x")
z <- Variable(1, name = "z")
obj <- Minimize(x + z)
copy_obj <- expr_copy(obj)
