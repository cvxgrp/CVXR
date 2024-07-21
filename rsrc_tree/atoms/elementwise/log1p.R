## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/log1p.py
#'
#' The Log1p class.
#'
#' This class represents the elementwise operation \eqn{\log(1 + x)}.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name Log1p-class
#' @aliases Log1p
#' @rdname Log1p-class
.Log1p <- setClass("Log1p", contains = "Log")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname Log1p-class
Log1p <- function(x) { .Log1p(x = x) }

#' @param object A \linkS4class{Log1p} object.
#' @param values A list of arguments to the atom.
#' @describeIn Log1p The elementwise natural logarithm of one plus the input value.
setMethod("to_numeric", "Log1p", function(object, values) { log(1 + values[[1]]) })

#' @describeIn Log1p The sign of the atom.
setMethod("sign_from_args", "Log1p", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @param values A list of numeric values for the arguments
#' @describeIn Log1p Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Log1p", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  # Outside domain or on boundary.
  if(min(values[[1]]) <= -1)
    return(list(NA_real_))   # Non-differentiable.
  else {
    grad_vals <- 1.0/(values[[1]] + 1)
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

#' @describeIn Log1p Returns constraints describng the domain of the node
setMethod(".domain", "Log1p", function(object) { list(object@args[[1]] >= -1) })

