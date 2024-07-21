## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/scalene.py

#'
#' The Scalene atom.
#'
#' This is an alias for \eqn{\alpha*Pos(x) + \beta*Neg(x)}.
#'
#' @param x An \linkS4class{Expression} object.
#' @param alpha A numeric constant.
#' @param beta A numeric constant.
#' @return An expression representing the scalene.
Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }

