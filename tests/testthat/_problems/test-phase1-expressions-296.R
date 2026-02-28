# Extracted from test-phase1-expressions.R:296

# test -------------------------------------------------------------------------
sm <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 3), x = c(1.0, 2.0),
                             dims = c(3, 3))
c1 <- Constant(sm)
canon <- canonicalize(c1)
