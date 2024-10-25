## CVXPY SOURCE: cvxpy/atoms/gen_lambda_max.py

#'
#' The GenLambdaMax class.
#'
#' This class represents the maximum generalized eigenvalue  \eqn{\lambda_{\max}(A, B)},
#' where \eqn{A} is a symmetric matrix and \eqn{B} is a positive semidefinite matrix.
#'
#' @slot A An \linkS4class{Expression} representing a symmetric matrix.
#' @slot B An \linkS4class{Expression} representing a positive semidefinite matrix.
#' @name GenLambdaMax-class
#' @aliases GenLambdaMax
#' @rdname GenLambdaMax-class
.GenLambdaMax <- setClass("GenLambdaMax", representation(A = "ConstValORExpr", B = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param B An \linkS4class{Expression} or numeric matrix.
#' @rdname GenLambdaMax-class
GenLambdaMax <- function(A, B) { .GenLambdaMax(A = A, B = B) }

setMethod("initialize", "GenLambdaMax", function(.Object, ..., A, B) {
  .Object@A <- A
  .Object@B <- B
  callNextMethod(.Object, ..., atom_args = list(.Object@A, .Object@B))
})

#' @param object A \linkS4class{GenLambdaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn GenLambdaMax The largest generalized eigenvalue corresponding to the matrices.
setMethod("to_numeric", "GenLambdaMax", function(object, values) {
  eigen_res <- geigen(values[[1]], values[[2]], symmetric = TRUE, only.values = TRUE)
  base::max(eigen_res$values)
})

#' @describeIn GenLambdaMax Returns constraints describing the domain of the node.
setMethod(".domain", "GenLambdaMax", function(object) { list(Conj(t(object@args[[1]])) == object@args[[1]], Conj(t(object@args[[2]])) == object@args[[2]], object@args[[2]] %>>% 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn GenLambdaMax Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "GenLambdaMax", function(object, values) { stop("Unimplemented") })

#' @describeIn GenLambdaMax Check that the matrices are square and of the same dimension.
setMethod("validate_args", "GenLambdaMax", function(object) {
  A_dim <- dim(object@args[[1]])
  B_dim <- dim(object@args[[2]])
  if(length(A_dim) != 2 || A_dim[1] != A_dim[2] || B_dim[1] != B_dim[2] || !all(A_dim == B_dim))
    stop("The arguments to GenLambdaMax must be square and have the same dimensions")
})

#' @describeIn GenLambdaMax The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "GenLambdaMax", function(object) { c(1,1) })

#' @describeIn GenLambdaMax The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "GenLambdaMax", function(object) { c(FALSE, FALSE) })

#' @describeIn GenLambdaMax Is the atom convex?
setMethod("is_atom_convex", "GenLambdaMax", function(object) { FALSE })

#' @describeIn GenLambdaMax Is the atom concave?
setMethod("is_atom_concave", "GenLambdaMax", function(object) { FALSE })

#' @describeIn GenLambdaMax Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "GenLambdaMax", function(object) { TRUE })

#' @describeIn GenLambdaMax Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "GenLambdaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn GenLambdaMax Is the atom weakly increasing in the index?
setMethod("is_incr", "GenLambdaMax", function(object, idx) { FALSE })

#' @describeIn GenLambdaMax Is the atom weakly decreasing in the index?
setMethod("is_decr", "GenLambdaMax", function(object, idx) { FALSE })

