## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/inv_pos.py

# x^{-1} for x > 0.
InvPos <- function(x) { Power(x, -1) }

