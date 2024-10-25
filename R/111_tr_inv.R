## CVXPY SOURCE: cvxpy/atoms/tr_inv.py
#'
#' The TrInv atom.
#'
#' The trace of the inverse of a positive definite matrix \eqn{X}, \eqn{tr(X^{-1})}.
#'
#' @slot X An \linkS4class{Expression} representing a positive definite matrix.
#' @name TrInv-class
#' @aliases TrInv
#' @rdname TrInv-class
.TrInv <- setClass("TrInv", representation(X = "ConstValORExpr"), contains = "Atom")

#' @param X An \linkS4class{Expression}
#' @rdname TrInv-class
TrInv <- function(X) { .TrInv(X = X) }

setMethod("initialize", "TrInv", function(.Object, ..., X) {
  .Object@X <- X
  callNextMethod(.Object, ..., atom_args = list(.Object@X))
})

#' @param object A \linkS4class{TrInv} object.
#' @param values A list of arguments to the atom.
#' @describeIn TrInv The trace of the inverse of a positive definite matrix.
setMethod("to_numeric", "TrInv", function(object, values) {
  # If values[[1]] isn't Hermitian, return Inf.
  if(base::norm(values[[1]] - Conj(t(values[[1]])), "F") >= 1e-8)
    return(Inf)
  # Take symmetric part of the input to enhance numerical stability.
  symm <- (values[[1]] + t(values[[1]]))/2
  eigVal <- eigen(symm, symmetric = TRUE, only.values = TRUE)$values
  if(min(eigVal) <= 0)
    return(Inf)
  return(base::sum(eigVal^-1))
})

#' @describeIn TrInv Check that the matrix is square.
setMethod("validate_args", "TrInv", function(object) {
  X <- object@args[[1]]
  if(len(dim(X)) == 1 || nrow(X) != ncol(X))
    stop("The argument to TrInv must be a square matrix")
})

#' @describeIn TrInv The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "TrInv", function(object) { c(1,1) })

#' @describeIn TrInv The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "TrInv", function(object) { c(TRUE, FALSE) })

#' @describeIn TrInv Is the atom convex?
setMethod("is_atom_convex", "TrInv", function(object) { TRUE })

#' @describeIn TrInv Is the atom concave?
setMethod("is_atom_concave", "TrInv", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn TrInv Is the atom weakly increasing in the index?
setMethod("is_incr", "TrInv", function(object, idx) { FALSE })

#' @describeIn TrInv Is the atom weakly decreasing in the index?
setMethod("is_decr", "TrInv", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn TrInv Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "TrInv", function(object, values) {
  X <- values[[1]]
  eigen_val <- eigen(X, only.values = TRUE)$values
  if(base::min(eigen_val) > 0) {
    # Grad: -t(X^{-2})
    D <- t(base::solve(X))
    D <- -D %*% D
    return(list(t(sparseMatrix(as.vector(D)))))
  } else   # Outside domain
    return(list(NA_real_))
})

#' @describeIn TrInv Returns constraints describing the domain of the node.
setMethod(".domain", "TrInv", function(object) { list(object@args[[1]] %>>% 0) })

#' @describeIn TrInv Returns the value of the atom.
setMethod("value", "TrInv", function(object) {
  arg_val <- value(object@args[[1]])
  if(base::abs(arg_val - Conj(t(arg_val))) > (ATOM_EVAL_TOL + ATOM_EVAL_TOL*base::abs(Conj(t(arg_val)))))
    stop("Input matrix was not Hermitian/symmetric")
  for(p in parameters(object)) {
    p_val <- value(p)
    if(is.na(p_val) || is.null(p_val))
      return(NA_real_)
  }
  return(value_impl(object))
})
