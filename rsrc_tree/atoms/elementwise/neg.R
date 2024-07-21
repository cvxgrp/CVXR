## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/neg.py
#'
#' An alias for -MinElemwise(x, 0)
#'
#' @param x An R numeric value or \linkS4class{Expression}.
#' @return An alias for -MinElemwise(x, 0)
#' @rdname Neg-int
Neg <- function(x) { -MinElemwise(x, 0) }
