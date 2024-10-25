## CVXPY SOURCE: cvxpy/atoms/suppfunc.py

#'
#' The SuppFuncAtom class.
#'
#' This class represents the support function of an affine expression \eqn{y}.
#'
#' @slot y An \linkS4class{Expression} with \code{is_affine(y) == TRUE}.
#' @slot parent A \linkS4class{SuppFunc} containing data for the convex set associated with this atom.
#' @name SuppFuncAtom-class
#' @aliases SuppFuncAtom
#' @rdname SuppFuncAtom-class
.SuppFuncAtom <- setClass("SuppFuncAtom", representation(y = "ConstValORExpr", parent = "SuppFunc", .eta = "numeric", .dim = "numeric"), contains = "Atom")

#' @param y An \linkS4class{Expression} with \code{is_affine(y) == TRUE}.
#' @param parent A \linkS4class{SuppFunc} containing data for the convex set associated with this atom.
#' @rdname SuppFuncAtom-class
SuppFuncAtom <- function(y, parent) { .SuppFuncAtom(y = y, parent = parent) }

setMethod("initialize", "SuppFuncAtom", function(.Object, ..., y, parent) {
  id(.Object) <- get_id()
  # .Object@args <- list(as.Constant(y))
  .Object@y <- y
  .Object@parent <- parent
  .Object@.eta <- NA_real_   # Store for debugging purposes.
  .Object@.dim <- c()
  validate_args()
  callNextMethod(.Object, ..., atom_args = list(as.Constant(y)))
})

#' @describeIn SuppFuncAtom List of \linkS4class{Variable} objects in the atom.
setMethod("variables", "SuppFuncAtom", function(object) { variables(object@args[[1]]) })

#' @describeIn SuppFuncAtom List of \linkS4class{Parameter} objects in the atom.
setMethod("parameters", "SuppFuncAtom", function(object) { list() })

#' @describeIn SuppFuncAtom List of \linkS4class{Constant} objects in the atom.
setMethod("constants", "SuppFuncAtom", function(object) { list() })

#' @param object A \linkS4class{SuppFuncAtom} object.
#' @describeIn SuppFuncAtom The dimensions of the atom.
setMethod("dim_from_args", "SuppFuncAtom", function(object) { object@.dim })

#' @describeIn SuppFuncAtom The sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "SuppFuncAtom", function(object) { c(FALSE, FALSE) })

#' @describeIn SuppFuncAtom A logical value indicating whether the atom is nonnegative.
setMethod("is_nonneg", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom A logical value indicating whether the atom is nonpositive.
setMethod("is_nonpos", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom A logical value indicating whether the atom is imaginary.
setMethod("is_imag", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom A logical value indicating whether the atom is complex valued.
setMethod("is_complex", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom convex?
setMethod("is_atom_convex", "SuppFuncAtom", function(object) { TRUE })

#' @describeIn SuppFuncAtom Is the atom concave?
setMethod("is_atom_concave", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "SuppFuncAtom", function(object) { TRUE })

#' @describeIn SuppFuncAtom Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "SuppFuncAtom", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SuppFuncAtom Is the atom weakly increasing in argument \code{idx}?
setMethod("is_incr", "SuppFuncAtom", function(object, idx) { FALSE })

#' @describeIn SuppFuncAtom Is the atom weakly decreasing in argument \code{idx}?
setMethod("is_decr", "SuppFuncAtom", function(object, idx) { FALSE })

#' @describeIn SuppFuncAtom Is the atom convex?
setMethod("is_convex", "SuppFuncAtom", function(object) { TRUE })

#' @describeIn SuppFuncAtom Is the atom concave?
setMethod("is_concave", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom log-log convex?
setMethod("is_log_log_convex", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom log-log concave?
setMethod("is_log_log_concave", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Is the atom quasiconvex?
setMethod("is_quasiconvex", "SuppFuncAtom", function(object) { TRUE })

#' @describeIn SuppFuncAtom Is the atom quasiconcave?
setMethod("is_quasiconcave", "SuppFuncAtom", function(object) { FALSE })

#' @describeIn SuppFuncAtom Returns the value of each of the components in the atom.
setMethod("value_impl", "SuppFuncAtom", function(object) {
  y_val <- value(object@args[[1]])
  y_val <- round(y_val, digits = 9)
  y_val <- as.vector(y_val)
  x_flat <- as.vector(t(object@parent@x))
  cons <- constraints(object@parent)
  if(len(cons) == 0) {
    dummy <- Variable()
    cons <- list(dummy == 1)
  }

  prob <- Problem(Maximize(y_val %*% x_flat), cons)
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  return(result$value)
})

#' @param values A list of numeric values for the arguments
#' @describeIn SuppFuncAtom Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SuppFuncAtom", function(object, values) {
  # The implementation of .grad from LogDet was used as a reference for this implementation.
  y0 <- value(object@args[[1]])   # Save for later.
  y <- values[[1]]   # Ignore all other values.
  value(object@args[[1]]) <- y
  .value_impl(object)   # TODO: What does dead-store do? Check Python functionality in suppfunc.py.
  value(object@args[[1]]) <- y0
  gradval <- value(object@parent@x)
  if(any(is.na(gradval)))
    # If we evaluated the support function successfully, then this means the support function is not finite at this point.
    return(list(NA_real_))
  else {
    gradmat <- t(sparseMatrix(as.vector(gradval)))
    return(gradmat)
  }
})

#' @describeIn SuppFuncAtom Strict inequalities are not allowed.
setMethod("<",  signature(e1 = "SuppFuncAtom", e2 = "Expression"), function(e1, e2) { stop("Unimplemented: Strict inequalities are not allowed.") })

#' @describeIn SuppFuncAtom Strict inequalities are not allowed.
setMethod(">",  signature(e1 = "SuppFuncAtom", e2 = "Expression"), function(e1, e2) { stop("Unimplemented: Strict inequalities are not allowed.") })
