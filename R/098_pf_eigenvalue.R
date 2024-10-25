## CVXPY SOURCE: cvxpy/atoms/pf_eigenvalue.py
#'
#' The PfEigenvalue class.
#'
#' This class represents the Perron-Frobenius eigenvalue of a positive matrix.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @name PfEigenvalue-class
#' @aliases PfEigenvalue
#' @rdname PfEigenvalue-class
.PfEigenvalue <- setClass("PfEigenvalue", representation(X = "ConstValORExpr"),
                          validity = function(object) {
                            if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                              stop("[PfEigenvalue: X] The argument X must be a square matrix")
                            return(TRUE)
                          }, contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname PfEigenvalue-class
PfEigenvalue <- function(X) { .PfEigenvalue(X = X) }

setMethod("initialize", "PfEigenvalue", function(.Object, ..., X = X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param x,object A \linkS4class{PfEigenvalue} object.
#' @describeIn PfEigenvalue The name and arguments of the atom.
setMethod("name", "PfEigenvalue", function(x) { paste(class(x), x@args[[1]]) })

#' @param values A list of arguments to the atom.
#' @describeIn PfEigenvalue Returns the Perron-Frobenius eigenvalue of \code{X}.
setMethod("to_numeric", "PfEigenvalue", function(object, values) {
  eig <- eigen(values[[1]], only.values = TRUE)
  max(abs(eig$values))
})

#' @describeIn PfEigenvalue The dimensions of the atom.
setMethod("dim_from_args", "PfEigenvalue", function(object) { c(1,1) })

#' @describeIn PfEigenvalue Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "PfEigenvalue", function(object) { c(TRUE, FALSE) })

#' @describeIn PfEigenvalue Is the atom convex?
setMethod("is_atom_convex", "PfEigenvalue", function(object) { FALSE })

#' @describeIn PfEigenvalue Is the atom concave?
setMethod("is_atom_concave", "PfEigenvalue", function(object) { FALSE })

#' @describeIn PfEigenvalue Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "PfEigenvalue", function(object) { TRUE })

#' @describeIn PfEigenvalue Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "PfEigenvalue", function(object) { FALSE })

# TODO: Figure out monotonicity.
#' @param idx An index into the atom.
#' @describeIn PfEigenvalue Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "PfEigenvalue", function(object, idx) { FALSE })

#' @param idx An index into the atom.
#' @describeIn PfEigenvalue Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "PfEigenvalue", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn PfEigenvalue Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "PfEigenvalue", function(object, values) { NA_real_ })
