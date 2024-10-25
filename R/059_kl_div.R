## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/kl_div.py
#'
#' The KLDiv class.
#'
#' The elementwise KL-divergence \eqn{x\log(x/y) - x + y}.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @slot y An \linkS4class{Expression} or numeric constant.
#' @name KLDiv-class
#' @aliases KLDiv
#' @rdname KLDiv-class
.KLDiv <- setClass("KLDiv", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @param y An \linkS4class{Expression} or numeric constant.
#' @rdname KLDiv-class
KLDiv <- function(x, y) { .KLDiv(x = x, y = y) }

setMethod("initialize", "KLDiv", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@y))
})

#' @param object A \linkS4class{KLDiv} object.
#' @param values A list of arguments to the atom.
#' @describeIn KLDiv The KL-divergence evaluted elementwise on the input value.
setMethod("to_numeric", "KLDiv", function(object, values) {
  x <- intf_convert_if_scalar(values[[1]])
  y <- intf_convert_if_scalar(values[[2]])

  # TODO: Return Inf outside domain
  xlogy <- function(x, y) {
    tmp <- x*log(y)
    tmp[x == 0] <- 0
    tmp
  }
  xlogy(x, x/y) - x + y
})

#' @describeIn KLDiv The atom is positive.
setMethod("sign_from_args", "KLDiv", function(object) { c(TRUE, FALSE) })

#' @describeIn KLDiv The atom is convex.
setMethod("is_atom_convex", "KLDiv", function(object) { TRUE })

#' @describeIn KLDiv The atom is not concave.
setMethod("is_atom_concave", "KLDiv", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn KLDiv The atom is not monotonic in any argument.
setMethod("is_incr", "KLDiv", function(object, idx) { FALSE })

#' @describeIn KLDiv The atom is not monotonic in any argument.
setMethod("is_decr", "KLDiv", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn KLDiv Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "KLDiv", function(object, values) {
  if(min(values[[1]]) <= 0 || min(values[[2]]) <= 0)
    return(list(NA_real_, NA_real_))   # Non-differentiable
  else {
    div <- values[[1]]/values[[2]]
    grad_vals <- list(log(div), 1-div)
    grad_list <- list()
    for(idx in 1:length(values)) {
      rows <- size(object@args[[idx]])
      cols <- size(object)
      grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals[[idx]], rows, cols)))
    }
    return(grad_list)
  }
})

#' @describeIn KLDiv Returns constraints describng the domain of the node
setMethod(".domain", "KLDiv", function(object) { list(object@args[[1]] >= 0, object@args[[2]] >= 0) })
