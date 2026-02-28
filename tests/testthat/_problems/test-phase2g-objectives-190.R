# Extracted from test-phase2g-objectives.R:190

# test -------------------------------------------------------------------------
x <- Variable(1, name = "x")
z <- Variable(1, name = "z")
obj <- Maximize(x + z)
copy_obj <- expr_copy(obj)
