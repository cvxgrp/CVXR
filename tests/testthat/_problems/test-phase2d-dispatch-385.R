# Extracted from test-phase2d-dispatch.R:385

# test -------------------------------------------------------------------------
skip_if_not_installed("Matrix")
x <- Variable(3)
A <- Matrix::Matrix(1:3, 3, 1)
result <- as.matrix(A) + x
expect_true(S7_inherits(result, AddExpression))
