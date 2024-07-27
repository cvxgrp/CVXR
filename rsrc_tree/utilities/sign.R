## CVXPY SOURCE: cvxpy/utilities/sign.py
###############################
#                             #
# Utility functions for signs #
#                             #
###############################
sum_signs <- function(exprs) {
  # Give the sign resulting from summing a list of expressions.
  is_pos <- all(sapply(exprs, is_nonneg))
  is_neg <- all(sapply(exprs, is_nonpos))
  c(is_pos, is_neg)
}

mul_sign <- function(lh_expr, rh_expr) {
  # Give the sign resulting from multiplying two expressions.
  # ZERO * ANYTHING == ZERO
  # POSITIVE * POSITIVE == POSITIVE
  # NEGATIVE * POSITIVE == NEGATIVE
  # NEGATIVE * NEGATIVE == POSITIVE

  lh_nonneg <- is_nonneg(lh_expr)
  rh_nonneg <- is_nonneg(rh_expr)
  lh_nonpos <- is_nonpos(lh_expr)
  rh_nonpos <- is_nonpos(rh_expr)

  lh_zero <- lh_nonneg && lh_nonpos
  rh_zero <- rh_nonneg && rh_nonpos

  is_zero <- lh_zero || rh_zero

  is_pos <- is_zero || (lh_nonneg && rh_nonneg) || (lh_nonpos && rh_nonpos)
  is_neg <- is_zero || (lh_nonneg && rh_nonpos) || (lh_nonpos && rh_nonneg)
  c(is_pos, is_neg)
}


