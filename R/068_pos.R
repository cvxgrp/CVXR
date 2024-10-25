## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/pos.py

#'
#' An alias for MaxElemwise(x, 0)
#'
#' @param x An R numeric value or \linkS4class{Expression}.
#' @return An alias for MaxElemwise(x, 0)
#' @rdname Pos-int
Pos <- function(x) { MaxElemwise(x, 0) }

