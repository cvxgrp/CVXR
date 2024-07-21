## CVXPY SOURCE: cvxpy/atoms/sum_smallest.py
#'
#' The SumSmallest atom.
#'
#' The sum of the smallest k values of a matrix.
#'
#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param k The number of smallest values to sum over.
#' @return Sum of the smlalest k values
SumSmallest <- function(x, k) {
  x <- as.Constant(x)
  -SumLargest(x = -x, k = k)
}

