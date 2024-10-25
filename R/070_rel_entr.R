## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/rel_entr.py
#'
#' The RelEntr class.
#'
#' This class represents elementwise \eqn{x\log(x/y)}.
#'
#' For disambiguation between RelEntr and KLDiv, see https://github.com/cvxpy/cvxpy/issues/733
#'
#' @slot x An \linkS4class{Expression}.
#' @slot y An \linkS4class{Expression}.
#' @name RelEntr-class
#' @aliases RelEntr
#' @rdname RelEntr-class
.RelEntr <- setClass("RelEntr", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression}.
#' @param y An \linkS4class{Expression}.
#' @rdname RelEntr-class
RelEntr <- function(x, y) { .RelEntr(x = x, y = y) }

setMethod("initialize", "RelEntr", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@y))
})

#' @param object A \linkS4class{RelEntr} object.
#' @param values A list of arguments to the atom.
#' @describeIn RelEntr The numeric value of the expression.
setMethod("to_numeric", "RelEntr", function(object, values) {
  x <- values[[1]]
  y <- values[[2]]
  rel_entr <- function(x_i, y_i) {
    if(x_i > 0 && y_i > 0)
      return(x_i*base::log(x_i/y_i))
    else if(x_i == 0 && y_i >= 0)
      return(0)
    else
      return(Inf)
  }
  return(mapply(rel_entr, x, y))
})

#' @describeIn RelEntr The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "RelEntr", function(object) { c(FALSE, FALSE) })

#' @describeIn RelEntr Is the atom convex?
setMethod("is_atom_convex", "RelEntr", function(object) { TRUE })

#' @describeIn RelEntr Is the atom concave?
setMethod("is_atom_concave", "RelEntr", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn RelEntr Is the atom weakly increasing in the index?
setMethod("is_incr", "RelEntr", function(object, idx) { FALSE })

#' @describeIn RelEntr Is the atom weakly decreasing in the index?
setMethod("is_decr", "RelEntr", function(object, idx) {
  if(idx == 1)
    return(FALSE)
  else
    return(TRUE)
})

#' @param values A list of numeric values for the arguments
#' @describeIn RelEntr Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "RelEntr", function(object, values) {
  if(min(values[[1]]) <= 0 || min(values[[2]]) <= 0)
    # Non-differentiable.
    return(list(NA_real_, NA_real_))
  else {
    div <- values[[1]]/values[[2]]
    grad_vals <- list(base::log(div) + 1, -div)
    grad_list <- list()
    for(idx in seq_along(values)) {
      rows <- size(object@args[[idx]])
      cols <- size(object)
      grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals[[idx]], rows, cols)))
    }
    return(grad_list)
  }
})

#' @describeIn RelEntr Returns constraints describing the domain of the node
setMethod(".domain", "RelEntr", function(object) { list(object@args[[1]] >= 0, object@args[[2]] >= 0) })
