## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/loggamma.py
#'
#' The LogGamma atom.
#'
#' The elementwise log of the gamma function.
#'
#' This implementation has modest accuracy over the full range, approaching perfect
#' accuracy as \eqn{x} approaches infinity. For details on the nature of the approximation,
#' return to CVXPY GitHub Issue #228 https://github.com/cvxpy/cvxpy/issues/228#issuecomment-544281906.
#'
#' @param x An \linkS4class{Expression} object.
#' @return An expression representing the log-gamma of \code{x}.
LogGamma <- function(x) {
  MaxElemwise(2.18382 - 3.62887*x, 1.79241 - 2.4902*x, 1.21628 - 1.37035*x,
              0.261474 - 0.28904*x, 0.577216 - 0.577216*x, -0.175517 + 0.03649*x,
              -1.27572 + 0.621514*x, -0.845568 + 0.422784*x, -0.577216*x - Log(x),
              0.918939 - x - Entr(x) - 0.5*Log(x))
}

