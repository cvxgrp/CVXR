## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/minimum.py

#'
#' The MinElemwise class.
#'
#' This class represents the elementwise minimum.
#'
#' @slot arg1 The first \linkS4class{Expression} in the minimum operation.
#' @slot arg2 The second \linkS4class{Expression} in the minimum operation.
#' @slot ... Additional \linkS4class{Expression} objects in the minimum operation.
#' @name MinElemwise-class
#' @aliases MinElemwise
#' @rdname MinElemwise-class
.MinElemwise <- setClass("MinElemwise", validity = function(object) {
                  if(is.null(object@args) || length(object@args) < 2)
                    stop("[MinElemwise: validation] args must have at least 2 arguments")
                  return(TRUE)
                }, contains = "Elementwise")

#' @param arg1 The first \linkS4class{Expression} in the minimum operation.
#' @param arg2 The second \linkS4class{Expression} in the minimum operation.
#' @param ... Additional \linkS4class{Expression} objects in the minimum operation.
#' @rdname MinElemwise-class
MinElemwise <- function(arg1, arg2, ...) { .MinElemwise(atom_args = list(arg1, arg2, ...)) }

#' @param object A \linkS4class{MinElemwise} object.
#' @param values A list of arguments to the atom.
#' @describeIn MinElemwise The elementwise minimum.
setMethod("to_numeric", "MinElemwise", function(object, values) {
  # Reduce(function(x, y) { ifelse(x <= y, x, y) }, values)
  Reduce("pmin", values)
})

#' @describeIn MinElemwise The sign of the atom.
setMethod("sign_from_args", "MinElemwise", function(object) {
  is_pos <- all(sapply(object@args, is_nonneg))
  is_neg <- any(sapply(object@args, is_nonpos))
  c(is_pos, is_neg)
})

#' @describeIn MinElemwise The atom is not convex.
setMethod("is_atom_convex", "MinElemwise", function(object) { FALSE })

#' @describeIn MinElemwise The atom is not concave.
setMethod("is_atom_concave", "MinElemwise", function(object) { TRUE })

#' @describeIn MinElemwise Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MinElemwise", function(object) { FALSE })

#' @describeIn MinElemwise Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MinElemwise", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinElemwise The atom is weakly increasing.
setMethod("is_incr", "MinElemwise", function(object, idx) { TRUE })

#' @describeIn MinElemwise The atom is not weakly decreasing.
setMethod("is_decr", "MinElemwise", function(object, idx) { FALSE })

#' @describeIn MinElemwise Are all the arguments piecewise linear?
setMethod("is_pwl", "MinElemwise", function(object) { all(sapply(object@args, is_pwl)) })

#' @param values A list of numeric values for the arguments
#' @describeIn MinElemwise Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MinElemwise", function(object, values) {
  min_vals <- to_numeric(object, values)
  vals_dim <- dim(min_vals)
  if(is.null(vals_dim))
    unused <- matrix(TRUE, nrow = length(min_vals), ncol = 1)
  else
    unused <- array(TRUE, dim = vals_dim)
  grad_list <- list()
  idx <- 1
  for(value in values) {
    rows <- size(object@args[[idx]])
    cols <- size(object)
    grad_vals <- (value == min_vals) & unused

    # Remove all the min_vals that were used
    unused[value == min_vals] <- FALSE
    grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
    idx <- idx + 1
  }
  grad_list
})
