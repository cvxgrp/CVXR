## CVXPY SOURCE: cvxpy/atoms/sigma_max.py
#'
#' The SigmaMax class.
#'
#' The maximum singular value of a matrix.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name SigmaMax-class
#' @aliases SigmaMax
#' @rdname SigmaMax-class
.SigmaMax <- setClass("SigmaMax", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or matrix.
#' @rdname SigmaMax-class
SigmaMax <- function(A = A) { .SigmaMax(A = A) }

setMethod("initialize", "SigmaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{SigmaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn SigmaMax The largest singular value of \code{A}.
setMethod("to_numeric", "SigmaMax", function(object, values) { base::norm(values[[1]], type = "2") })

#' @describeIn SigmaMax Does the atom handle complex numbers?
setMethod("allow_complex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is a scalar.
setMethod("dim_from_args", "SigmaMax", function(object) { c(1,1) })

#' @describeIn SigmaMax The atom is positive.
setMethod("sign_from_args",  "SigmaMax", function(object) { c(TRUE, FALSE) })

#' @describeIn SigmaMax The atom is convex.
setMethod("is_atom_convex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is concave.
setMethod("is_atom_concave", "SigmaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_incr", "SigmaMax", function(object, idx) { FALSE })

#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_decr", "SigmaMax", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn SigmaMax Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SigmaMax", function(object, values) {
  # Grad: U diag(e_1) t(V)
  s <- svd(values[[1]])
  ds <- rep(0, length(s$d))
  ds[1] <- 1
  D <- s$u %*% diag(ds) %*% t(s$v)
  list(Matrix(as.vector(D), sparse = TRUE))
})

