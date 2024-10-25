## CVXPY SOURCE: cvxpy/atoms/affine/vec.py

# Reshape into single column vector.
Vec <- function(X, byrow = FALSE) {
  X <- as.Constant(X)
  # Reshape(expr = X, new_dim = size(X), byrow = byrow)
  Reshape(expr = X, new_dim = c(size(X), 1), byrow = byrow)
}
