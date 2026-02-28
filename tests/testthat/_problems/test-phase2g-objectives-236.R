# Extracted from test-phase2g-objectives.R:236

# test -------------------------------------------------------------------------
m <- Matrix::sparseMatrix(i = 1, j = 1, x = 3, dims = c(1, 1))
expect_equal(scalar_value(m), 3)
