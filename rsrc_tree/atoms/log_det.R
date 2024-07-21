## CVXPY SOURCE: cvxpy/atoms/log_det.py

#'
#' The LogDet class.
#'
#' The natural logarithm of the determinant of a matrix, \eqn{\log\det(A)}.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name LogDet-class
#' @aliases LogDet
#' @rdname LogDet-class
.LogDet <- setClass("LogDet", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname LogDet-class
LogDet <- function(A) { .LogDet(A = A) }

setMethod("initialize", "LogDet", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{LogDet} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogDet The log-determinant of SDP matrix \code{A}. This is the sum of logs of the eigenvalues and is equivalent to the nuclear norm of the matrix logarithm of \code{A}.
setMethod("to_numeric", "LogDet", function(object, values) {
  if(is.complex(values[[1]])) {
    eigvals <- eigen(values[[1]], only.values = TRUE)$values
    return(log(prod(eigvals)))
  } else {
    logdet <- determinant(values[[1]], logarithm = TRUE)
    if(logdet$sign == 1)
      return(as.numeric(logdet$modulus))
    else
      return(-Inf)
  }
})

#' @describeIn LogDet Check that \code{A} is square.
setMethod("validate_args", "LogDet", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(length(arg_dim) == 1 || arg_dim[1] != arg_dim[2])
    stop("The argument to LogDet must be a square matrix")
})

#' @describeIn LogDet The atom is a scalar.
setMethod("dim_from_args", "LogDet", function(object) { c(1,1) })

#' @describeIn LogDet The atom is non-negative.
setMethod("sign_from_args",  "LogDet", function(object) { c(TRUE, FALSE) })

#' @describeIn LogDet The atom is not convex.
setMethod("is_atom_convex", "LogDet", function(object) { FALSE })

#' @describeIn LogDet The atom is concave.
setMethod("is_atom_concave", "LogDet", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_incr", "LogDet", function(object, idx) { FALSE })

#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_decr", "LogDet", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn LogDet Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LogDet", function(object, values) {
  X <- as.matrix(values[[1]])
  eigen_val <- eigen(X, only.values = TRUE)$values
  if(min(eigen_val) > 0) {
    # Grad: t(X^(-1))
    D <- t(base::solve(X))
    return(list(Matrix(as.vector(D), sparse = TRUE)))
  } else   # Outside domain
    return(list(NA_real_))
})

#' @describeIn LogDet Returns constraints describing the domain of the node
setMethod(".domain", "LogDet", function(object) { list(object@args[[1]] %>>% 0) })

