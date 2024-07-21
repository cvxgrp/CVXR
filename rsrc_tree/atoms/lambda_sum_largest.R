## CVXPY SOURCE: cvxpy/atoms/lambda_sum_largest.py
#'
#' The LambdaSumLargest class.
#'
#' This class represents the sum of the \code{k} largest eigenvalues of a matrix.
#'
#' @slot k A positive integer.
#' @name LambdaSumLargest-class
#' @aliases LambdaSumLargest
#' @rdname LambdaSumLargest-class
.LambdaSumLargest <- setClass("LambdaSumLargest", representation(k = "numeric"), contains = "LambdaMax")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param k A positive integer.
#' @rdname LambdaSumLargest-class
LambdaSumLargest <- function(A, k) { .LambdaSumLargest(A = A, k = k) }

setMethod("initialize", "LambdaSumLargest", function(.Object, ..., k) {
  .Object@k <- k
  callNextMethod(.Object, ...)
})

#' @describeIn LambdaSumLargest Does the atom handle complex numbers?
setMethod("allow_complex", "LambdaSumLargest", function(object) { TRUE })

#' @param object A \linkS4class{LambdaSumLargest} object.
#' @param values A list of arguments to the atom.
#' @describeIn LambdaSumLargest Returns the largest eigenvalue of \code{A}, which must be symmetric.
setMethod("to_numeric", "LambdaSumLargest", function(object, values) {
  # if(any(t(values[[1]]) != values[[1]]))
  #  stop("LambdaSumLargest called on a non-symmetric matrix")
  eigs <- eigen(values[[1]], only.values = TRUE)$values
  value(SumLargest(eigs, object@k))
})

#' @describeIn LambdaSumLargest Verify that the argument \code{A} is square.
setMethod("validate_args", "LambdaSumLargest", function(object) {
  A <- object@args[[1]]
  if(ndim(A) != 2 || nrow(A) != ncol(A))
    stop("First argument must be a square matrix.")
  else if(as.integer(object@k) != object@k || object@k <= 0)
    stop("Second argument must be a positive integer.")
})

#' @describeIn LambdaSumLargest Returns the parameter \code{k}.
setMethod("get_data", "LambdaSumLargest", function(object) { list(object@k) })

#' @param values A list of numeric values for the arguments
#' @describeIn LambdaSumLargest Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LambdaSumLargest", function(object, values) { stop("Unimplemented") })
