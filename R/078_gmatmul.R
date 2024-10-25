## CVXPY SOURCE: cvxpy/atoms/gmatmul.py
#'
#' The GmatMul atom.
#'
#' Geometric matrix multiplication: \eqn{A \mathbin{\diamond} X}.
#'
#' For \eqn{A \in \mathbf{R}^{m \times n}} and \eqn{X \in \mathbf{R}^{n \times p}_{++}}, this atom represents
#'
#' \deqn{\left[\begin{array}{ccc} \prod_{j=1}^n X_{j1}^{A_{1j}} & \cdots & \prod_{j=1}^n X_{pj}^{A_{1j}} \\ \vdots &  & \vdots \\ \prod_{j=1}^n X_{j1}^{A_{mj}} & \cdots & \prod_{j=1}^n X_{pj}^{A_{mj}} \end{array}\right]}
#'
#' This atom is log-log affine in \eqn{X}.
#' @slot A An \linkS4class{Expression} representing a constant matrix.
#' @slot X An \linkS4class{Expression} representing a positive matrix.
#' @name GmatMul-class
#' @aliases GmatMul
#' @rdname GmatMul-class
.GmatMul <- setClass("GmatMul", representation(A = "ConstValORExpr", X = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname GmatMul-class
GmatMul <- function(A, X) { .GmatMul(A = A, X = X) }

setMethod("initialize", "GmatMul", function(.Object, ..., A, X) {
  .Object@A <- as.Constant(A)
  .Object@X <- X
  callNextMethod(.Object, ..., atom_args = list(.Object@X))
})

#' @param object A \linkS4class{GmatMul} object.
#' @param values A list of arguments to the atom.
#' @describeIn GmatMul Geometric matrix multiplication.
setMethod("to_numeric", "GmatMul", function(object, values) {
  logX <- base::log(values[[1]])
  return(base::exp(value(object@A) %*% logX))
})

#' @describeIn GmatMul The name and arguments of the atom.
setMethod("name", "GmatMul", function(x) {
  sprintf("%s(%s, %s)", class(x), name(x@A), name(x@args[[1]]))
})

#' @describeIn GmatMul Check if the arguments are valid.
setMethod("validate_args", "GmatMul", function(object) {
  callNextMethod()
  if(!is_constant(object@A))
    stop("A must be constant")
  if(length(parameters(A)) != 0 && !is(object@A, Parameter))
    stop("A must be of class Constant or Parameter")
  if(!is_pos(object@args[[1]]))
    stop("X must be positive")
})

#' @describeIn GmatMul The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "GmatMul", function(object) { mul_shapes(dim(object@A), dim(object@args[[1]])) })

#' @describeIn GmatMul Returns the parameter \code{A}.
setMethod("get_data", "GmatMul", function(object) { list(object@A) })

#' @describeIn GmatMul The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "GmatMul", function(object) { c(TRUE, FALSE) })

#' @describeIn GmatMul Is the atom convex?
setMethod("is_atom_convex", "GmatMul", function(object) { FALSE })

#' @describeIn GmatMul Is the atom concave?
setMethod("is_atom_concave", "GmatMul", function(object) { FALSE })

#' @describeIn GmatMul List of \linkS4class{Parameter} objects in the atom.
setMethod("parameters", "GmatMul", function(object) {
  # The exponent matrix, which is not an argument, may be parametrized.
  c(parameters(object@args[[1]]), parameters(obje))
})

#' @describeIn GmatMul Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "GmatMul", function(object) {
  if(dpp_scope_active()) {
    # This branch applies curvature rules for DPP.
    #
    # Because a DPP scope is active, parameters will be
    # treated as affine (like variables, not constants) by curvature
    # analysis methods.
    #
    # A power X^A is log-log convex (actually, affine) as long as
    # at least one of X and P do not contain parameters.
    #
    # Note by construction (see A is either a Constant or a Parameter, ie,
    # either is(A, Constant) or is(A, Parameter)).
    X <- object@args[[1]]
    A <- object@A
    return(!(length(parameters(X)) != 0 && length(parameters(A)) != 0))
  } else
    return(TRUE)
})

#' @describeIn GmatMul Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "GmatMul", function(object) { is_log_log_convex(object) })

#' @param idx An index into the atom.
#' @describeIn GmatMul Is the atom weakly increasing in the index?
setMethod("is_incr", "GmatMul", function(object, idx) { is_nonneg(object@A) })

#' @describeIn GmatMul Is the atom weakly decreasing in the index?
setMethod("is_decr", "GmatMul", function(object, idx) { is_nonpos(object@A) })

#' @param values A list of numeric values for the arguments
#' @describeIn GmatMul Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "GmatMul", function(object, values) { return(NA_real_) })

