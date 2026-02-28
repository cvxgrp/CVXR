# Extracted from test-phase5a-canonicalizers.R:290

# test -------------------------------------------------------------------------
x <- Variable(shape = c(3L, 1L))
y <- Variable(shape = c(1L, 1L), nonneg = TRUE)
expr <- QuadOverLin(x, y)
