## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/log_normcdf.py
#'
#' The LogNormCdf atom.
#'
#' The elementwise log of the cumulative distribution function of a standard normal random variable.
#'
#' This implementation is a quadratic approximation with modest accuracy over [-4, 4].
#' For details on the nature of the approximation, refer to CVXPY GitHub PR #1224
#' https://github.com/cvxpy/cvxpy/pull/1224#issue-793221374.
#'
#' Note: Python SciPy's analog of LogNormCdf is called log_ndtr:
#' https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_ndtr.html
#'
#' @param x An \linkS4class{Expression} object.
#' @return An expression representing the log-norm of \code{x}.
LogNormCdf <- function(x) {
  x <- as.Constant(x)
  flat_x <- Reshape(x, c(1, size(x)))

  d <- sqrt(c(0.02301291, 0.08070214, 0.16411522, 0.09003495, 0.08200854,
            0.01371543, 0.04641081))
  A <- sparseMatrix(i = seq(length(d)), j = seq(length(d)), x = d)
  b <- matrix(c(3, 2, 1, 0,-1, -2.5, -3.5))

  y <- A %*% (b %*% matrix(1, dim(flat_x)) - matrix(1, dim(b)) %*% flat_x)
  out <- -SumEntries(MaxElemwise(y, 0)^2, axis = 2)

  return(Reshape(out, dim(x)))
}
