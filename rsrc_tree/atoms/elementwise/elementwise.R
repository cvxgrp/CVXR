## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/elementwise.py

#'
#' The Elementwise class.
#'
#' This virtual class represents an elementwise atom.
#'
#' @name Elementwise-class
#' @aliases Elementwise
#' @rdname Elementwise-class
Elementwise <- setClass("Elementwise", contains = c("VIRTUAL", "Atom"))

#' @param object An \linkS4class{Elementwise} object.
#' @describeIn Elementwise Dimensions is the same as the sum of the arguments' dimensions.
setMethod("dim_from_args", "Elementwise", function(object) {
  sum_dims(lapply(object@args, dim))
})

#' @describeIn Elementwise Verify that all the dimensions are the same or can be promoted.
setMethod("validate_args", "Elementwise", function(object) {
  sum_dims(lapply(object@args, dim))
  callNextMethod()
})

#' @describeIn Elementwise Is the expression symmetric?
setMethod("is_symmetric", "Elementwise", function(object) {
  symm_args <- all(sapply(object@args, is_symmetric))
  return(nrow(object) == ncol(object) && symm_args)
})

#
# Gradient to Diagonal
#
# Converts elementwise gradient into a diagonal matrix.
#
# @param value A scalar value or matrix.
# @return A sparse matrix.
# @rdname Elementwise-elemwise_grad_to_diag
Elementwise.elemwise_grad_to_diag <- function(value, rows, cols) {
  value <- as.vector(value)
  sparseMatrix(i = 1:rows, j = 1:cols, x = value, dims = c(rows, cols))
}

#
# Promotes LinOp
#
# Promotes the LinOp if necessary.
# @param arg The LinOp to promote.
# @param dim The desired dimensions.
# @return The promoted LinOp.
# @rdname Elementwise-promote
Elementwise.promote <- function(arg, dim) {
  if(any(dim(arg) != dim))
    lu.promote(arg, dim)
  else
    arg
}
