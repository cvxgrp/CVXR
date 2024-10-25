## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/power.py
#'
#' The Power class.
#'
#' This class represents the elementwise power function \eqn{f(x) = x^p}.
#' If \code{expr} is a CVXR expression, then \code{expr^p} is equivalent to \code{Power(expr, p)}.
#'
#' For \eqn{p = 0}, \eqn{f(x) = 1}, constant, positive.
#'
#' For \eqn{p = 1}, \eqn{f(x) = x}, affine, increasing, same sign as \eqn{x}.
#'
#' For \eqn{p = 2,4,8,...}, \eqn{f(x) = |x|^p}, convex, signed monotonicity, positive.
#'
#' For \eqn{p < 0} and \eqn{f(x) = }
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x > 0}}
#'   \item{\eqn{+\infty}}{\eqn{x \leq 0}}
#' }, this function is convex, decreasing, and positive.
#'
#' For \eqn{0 < p < 1} and \eqn{f(x) =}
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x \geq 0}}
#'   \item{\eqn{-\infty}}{\eqn{x < 0}}
#' }, this function is concave, increasing, and positive.
#'
#' For \eqn{p > 1, p \neq 2,4,8,\ldots} and \eqn{f(x) = }
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x \geq 0}}
#'   \item{\eqn{+\infty}}{\eqn{x < 0}}
#' }, this function is convex, increasing, and positive.
#'
#' @slot x The \linkS4class{Expression} to be raised to a power.
#' @slot p A numeric value indicating the scalar power.
#' @slot max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @name Power-class
#' @aliases Power
#' @rdname Power-class
.Power <- setClass("Power", representation(x = "ConstValORExpr", p = "NumORgmp", max_denom = "numeric", w = "NumORgmp", approx_error = "numeric"),
                          prototype(max_denom = 1024, w = NA_real_, approx_error = NA_real_), contains = "Elementwise")

#' @param x The \linkS4class{Expression} to be raised to a power.
#' @param p A numeric value indicating the scalar power.
#' @param max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @rdname Power-class
Power <- function(x, p, max_denom = 1024) { .Power(x = x, p = p, max_denom = max_denom) }

setMethod("initialize", "Power", function(.Object, ..., x, p, max_denom = 1024, w = NA_real_, approx_error = NA_real_) {
  p_old <- p
  if(length(p) != 1)
    stop("p must be a numeric scalar")
  if(is.na(p) || is.null(p))
    stop("p cannot be NA or NULL")

  # How we convert p to a rational depends on the branch of the function
  if(p > 1) {
    pw <- pow_high(p)
    p <- pw[[1]]
    w <- pw[[2]]
  } else if(p > 0 && p < 1) {
    pw <- pow_mid(p)
    p <- pw[[1]]
    w <- pw[[2]]
  } else if(p < 0) {
    pw <- pow_neg(p)
    p <- pw[[1]]
    w <- pw[[2]]
  }

  if(p == 1) {
    # In case p is a fraction equivalent to 1
    p <- 1
    w <- NA_real_
  } else if(p == 0) {
    p <- 0
    w <- NA_real_
  }

  .Object@p <- p
  .Object@w <- w
  .Object@approx_error <- as.double(abs(.Object@p - p_old))

  .Object@x <- x
  .Object@max_denom <- max_denom
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{Power} object.
#' @param values A list of arguments to the atom.
#' @describeIn Power Throw an error if the power is negative and cannot be handled.
setMethod("to_numeric", "Power", function(object, values) {
  # Throw error if negative and Power doesn't handle that
  if(object@p < 0 && min(values[[1]]) <= 0)
    stop("Power cannot be applied to negative or zero values")
  else if(!is_power2(object@p) && object@p != 0 && min(values[[1]]) < 0)
    stop("Power cannot be applied to negative values")
  else
    return(values[[1]]^(as.double(object@p)))
})

#' @describeIn Power The sign of the atom.
setMethod("sign_from_args", "Power", function(object) {
  if(object@p == 1)   # Same as input
    c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]]))
  else   # Always positive
    c(TRUE, FALSE)
})

#' @describeIn Power Is \eqn{p \leq 0} or \eqn{p \geq 1}?
setMethod("is_atom_convex", "Power", function(object) {
  # p == 0 is affine here.
  object@p <= 0 || object@p >= 1
})

#' @describeIn Power Is \eqn{p \geq 0} or \eqn{p \leq 1}?
setMethod("is_atom_concave", "Power", function(object) {
  # p == 0 is affine here.
  object@p >= 0 && object@p <= 1
})

#' @describeIn Power Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Power", function(object) { TRUE })

#' @describeIn Power Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Power", function(object) { TRUE })

#' @describeIn Power A logical value indicating whether the atom is constant.
setMethod("is_constant", "Power", function(object) { object@p == 0 || callNextMethod() })

#' @param idx An index into the atom.
#' @describeIn Power A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Power", function(object, idx) {
  if(object@p >= 0 && object@p <= 1)
    return(TRUE)
  else if(object@p > 1) {
    if(is_power2(object@p))
      return(is_nonneg(object@args[[idx]]))
    else
      return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Power A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Power", function(object, idx) {
  if(object@p <= 0)
    return(TRUE)
  else if(object@p > 1) {
    if(is_power2(object@p))
      return(is_nonpos(object@args[[idx]]))
    else
      return(FALSE)
  } else
    return(FALSE)
})

#' @describeIn Power A logical value indicating whether the atom is quadratic.
setMethod("is_quadratic", "Power", function(object) {
  if(object@p == 0)
    return(TRUE)
  else if(object@p == 1)
    return(is_quadratic(object@args[[1]]))
  else if(object@p == 2)
    return(is_affine(object@args[[1]]))
  else
    return(is_constant(object@args[[1]]))
})

#' @describeIn Power A logical value indicating whether the atom is quadratic of piecewise affine.
setMethod("is_qpwa", "Power", function(object) {
  if(object@p == 0)
    return(TRUE)
  else if(object@p == 1)
    return(is_qpwa(object@args[[1]]))
  else if(object@p == 2)
    return(is_pwl(object@args[[1]]))
  else
    return(is_constant(object@args[[1]]))
})

#' @param values A list of numeric values for the arguments
#' @describeIn Power Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Power", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  if(object@p == 0) # All zeros
    return(list(sparseMatrix(i = c(), j = c(), dims = c(rows, cols))))

  # Outside domain or on boundary
  if(!is_power2(object@p) && min(values[[1]]) <= 0) {
    if(object@p < 1)
      return(list(NA_real_))  # Non-differentiable
    else   # Round up to zero
      values[[1]] <- ifelse(values[[1]] >= 0, values[[1]], 0)
  }

  grad_vals <- as.double(object@p) * (values[[1]]^(as.double(object@p) - 1))
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

#' @describeIn Power Returns constraints describing the domain of the node
setMethod(".domain", "Power", function(object) {
  if((object@p < 1 && object@p != 0) || (object@p > 1 && !is_power2(object@p)))
    list(object@args[[1]] >= 0)
  else
    list()
})

#' @describeIn Power A list containing the output of \code{pow_low, pow_mid}, or \code{pow_high} depending on the input power.
setMethod("get_data", "Power", function(object) { list(object@p, object@w) })

#' @param args A list of arguments to reconstruct the atom. If args=NULL, use the current args of the atom
#' @param id_objects Currently unused.
#' @describeIn Power Returns a shallow  copy of the power atom
setMethod("copy", "Power", function(object, args = NULL, id_objects = list()) {
  if(is.null(args))
    args <- object@args
  data <- get_data(object)
  copy <- do.call(class(object), args)
  copy@p <- data[[1]]
  copy@w <- data[[2]]
  copy@approx_error <- object@approx_error
  copy
})

#' @describeIn Power Returns the expression in string form.
setMethod("name", "Power", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ", ", x@p, ")", sep = "")
})

