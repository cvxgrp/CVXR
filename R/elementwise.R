#'
#' The Elementwise class.
#'
#' This virtual class represents an elementwise atom.
#'
#' @name Elementwise-class
#' @aliases Elementwise
#' @rdname Elementwise-class
Elementwise <- setClass("Elementwise", contains = c("VIRTUAL", "Atom"))

#' @param object An \linkS4class{Elementwise} object.
#' @describeIn Elementwise Dimensions is the same as the sum of the arguments' dimensions.
setMethod("dim_from_args", "Elementwise", function(object) {
  sum_dims(lapply(object@args, dim))
})

#' @describeIn Elementwise Verify that all the dimensions are the same or can be promoted.
setMethod("validate_args", "Elementwise", function(object) {
  sum_dims(lapply(object@args, dim))
  callNextMethod()
})

#' @describeIn Elementwise Is the expression symmetric?
setMethod("is_symmetric", "Elementwise", function(object) {
  symm_args <- all(sapply(object@args, is_symmetric))
  return(nrow(object) == ncol(object) && symm_args)
})

#
# Gradient to Diagonal
#
# Converts elementwise gradient into a diagonal matrix.
#
# @param value A scalar value or matrix.
# @return A sparse matrix.
# @rdname Elementwise-elemwise_grad_to_diag
Elementwise.elemwise_grad_to_diag <- function(value, rows, cols) {
  value <- as.vector(value)
  sparseMatrix(i = 1:rows, j = 1:cols, x = value, dims = c(rows, cols))
}

#
# Promotes LinOp
#
# Promotes the LinOp if necessary.
# @param arg The LinOp to promote.
# @param dim The desired dimensions.
# @return The promoted LinOp.
# @rdname Elementwise-promote
Elementwise.promote <- function(arg, dim) {
  if(any(dim(arg) != dim))
    lo.promote(arg, dim)
  else
    arg
}

#'
#' The Abs class.
#'
#' This class represents the elementwise absolute value.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Abs-class
#' @aliases Abs
#' @rdname Abs-class
.Abs <- setClass("Abs", representation(x = "Expression"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} object.
#' @rdname Abs-class
Abs <- function(x) { .Abs(x = x) }

setMethod("initialize", "Abs", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object An \linkS4class{Abs} object.
#' @param values A list of arguments to the atom.
#' @describeIn Abs The elementwise absolute value of the input.
setMethod("to_numeric", "Abs", function(object, values) { abs(values[[1]]) })

#' @describeIn Abs Does the atom handle complex numbers?
setMethod("allow_complex", "Abs", function(object) { TRUE })

#' @describeIn Abs The atom is positive.
setMethod("sign_from_args", "Abs", function(object) { c(TRUE, FALSE) })

#' @describeIn Abs The atom is convex.
setMethod("is_atom_convex", "Abs", function(object) { TRUE })

#' @describeIn Abs The atom is not concave.
setMethod("is_atom_concave", "Abs", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Abs A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Abs", function(object, idx) { is_nonneg(object@args[[idx]]) })

#' @param idx An index into the atom.
#' @describeIn Abs A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Abs", function(object, idx) { is_nonpos(object@args[[idx]]) })

#' @describeIn Abs Is \code{x} piecewise linear?
setMethod("is_pwl", "Abs", function(object) { 
  is_pwl(object@args[[1]]) && (is_real(object@args[[1]]) || is_imag(object@args[[1]]))
})

setMethod(".grad", "Abs", function(object, values) {
  # Grad: +1 if positive, -1 if negative
  rows <- size(expr(object))
  cols <- size(object)
  D <- array(0, dim = dim(expr(object)))
  D <- D + (values[[1]] > 0)
  D <- D - (values[[1]] < 0)
  list(Elementwise.elemwise_grad_to_diag(D, rows, cols))
})

#'
#' The Entr class.
#'
#' This class represents the elementwise operation \eqn{-xlog(x)}.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name Entr-class
#' @aliases Entr
#' @rdname Entr-class
.Entr <- setClass("Entr", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname Entr-class
Entr <- function(x) { .Entr(x = x) }

setMethod("initialize", "Entr", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object An \linkS4class{Entr} object.
#' @param values A list of arguments to the atom.
#' @describeIn Entr The elementwise entropy function evaluated at the value.
setMethod("to_numeric", "Entr", function(object, values) {
  xlogy <- function(x, y) {
    tmp <- x*log(y)
    tmp[x == 0] <- 0
    tmp
  }

  x <- values[[1]]
  results <- -xlogy(x, x)

  # Return -Inf outside the domain
  results[is.na(results)] <- -Inf
  if(all(dim(results) == 1))
    results <- as.vector(results)
  results
})

#' @describeIn Entr The sign of the atom is unknown.
setMethod("sign_from_args", "Entr", function(object) { c(FALSE, FALSE) })

#' @describeIn Entr The atom is not convex.
setMethod("is_atom_convex", "Entr", function(object) { FALSE })

#' @describeIn Entr The atom is concave.
setMethod("is_atom_concave", "Entr", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn Entr The atom is weakly increasing.
setMethod("is_incr", "Entr", function(object, idx) { FALSE })

#' @describeIn Entr The atom is weakly decreasing.
setMethod("is_decr", "Entr", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Entr Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Entr", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  # Outside domain or on boundary
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- -log(values[[1]]) - 1
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

#' @describeIn Entr Returns constraints descrbing the domain of the node
setMethod(".domain", "Entr", function(object) { list(object@args[[1]] >= 0) })

#'
#' The Exp class.
#'
#' This class represents the elementwise natural exponential \eqn{e^x}.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Exp-class
#' @aliases Exp
#' @rdname Exp-class
.Exp <- setClass("Exp", representation(x = "Expression"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} object.
#' @rdname Exp-class
Exp <- function(x) { .Exp(x = x) }

setMethod("initialize", "Exp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object An \linkS4class{Exp} object.
#' @param values A list of arguments to the atom.
#' @describeIn Exp The matrix with each element exponentiated.
setMethod("to_numeric", "Exp", function(object, values) { exp(values[[1]]) })

#' @describeIn Exp The atom is positive.
setMethod("sign_from_args", "Exp", function(object) { c(TRUE, FALSE) })

#' @describeIn Exp The atom is convex.
setMethod("is_atom_convex", "Exp", function(object) { TRUE })

#' @describeIn Exp The atom is not concave.
setMethod("is_atom_concave", "Exp", function(object) { FALSE })

#' @describeIn Exp Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Exp", function(object) { TRUE })

#' @describeIn Exp Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Exp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Exp The atom is weakly increasing.
setMethod("is_incr", "Exp", function(object, idx) { TRUE })

#' @describeIn Exp The atom is not weakly decreasing.
setMethod("is_decr", "Exp", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Exp Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Exp", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)
  grad_vals <- exp(values[[1]])
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

#'
#' The Huber class.
#'
#' This class represents the elementwise Huber function, \eqn{Huber(x, M) = }
#' \itemize{
#'   \item{\eqn{2M|x|-M^2}}{for \eqn{|x| \geq |M|}}
#'    \item{\eqn{|x|^2}}{for \eqn{|x| \leq |M|.}}
#'  }
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @slot M A positive scalar value representing the threshold. Defaults to 1.
#' @name Huber-class
#' @aliases Huber
#' @rdname Huber-class
.Huber <- setClass("Huber", representation(x = "ConstValORExpr", M = "ConstValORExpr"),
                           prototype(M = 1), contains = "Elementwise")

#' @param x An \linkS4class{Expression} object.
#' @param M A positive scalar value representing the threshold. Defaults to 1.
#' @rdname Huber-class
Huber <- function(x, M = 1) { .Huber(x = x, M = M) }

setMethod("initialize", "Huber", function(.Object, ..., x, M = 1) {
  .Object@M <- as.Constant(M)
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{Huber} object.
#' @param values A list of arguments to the atom.
#' @describeIn Huber The Huber function evaluted elementwise on the input value.
setMethod("to_numeric", "Huber", function(object, values) {
  huber_loss <- function(delta, r) {
    if(delta < 0)
      return(Inf)
    else if(delta >= 0 && abs(r) <= delta)
      return(r^2/2)
    else
      return(delta * (abs(r) - delta/2))
  }

  M_val <- value(object@M)
  val <- values[[1]]
  if(is.null(dim(val)))
    result <- 2*huber_loss(M_val, val)
  else if(is.vector(val))
    result <- 2*sapply(val, function(v) { huber_loss(M_val, v) })
  else
    result <- 2*apply(val, 1:length(dim(val)), function(v) { huber_loss(M_val, v) })
  
  if(all(dim(result) == 1))
    result <- as.vector(result)
  return(result)
})

#' @describeIn Huber The atom is positive.
setMethod("sign_from_args", "Huber", function(object) { c(TRUE, FALSE) })

#' @describeIn Huber The atom is convex.
setMethod("is_atom_convex", "Huber", function(object) { TRUE })

#' @describeIn Huber The atom is not concave.
setMethod("is_atom_concave", "Huber", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Huber A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Huber", function(object, idx) { is_nonneg(object@args[[idx]]) })

#' @describeIn Huber A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Huber", function(object, idx) { is_nonpos(object@args[[idx]]) })

#' @describeIn Huber The atom is quadratic if \code{x} is affine.
setMethod("is_quadratic", "Huber", function(object) { is_affine(object@args[[1]]) })

#' @describeIn Huber A list containing the parameter \code{M}.
setMethod("get_data", "Huber", function(object) { list(object@M) })

#' @describeIn Huber Check that \code{M} is a non-negative constant.
setMethod("validate_args", "Huber", function(object) {
  if(!(is_nonneg(object@M) && is_constant(object@M) && is_scalar(object@M)))
    stop("M must be a non-negative scalar constant")
  callNextMethod()
})

#' @param values A list of numeric values for the arguments
#' @describeIn Huber Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Huber", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  val_abs <- abs(values[[1]])
  M_val <- as.numeric(value(object@M))
  min_val <- ifelse(val_abs >= M_val, M_val, val_abs)

  grad_vals <- 2*(sign(values[[1]]) * min_val)
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

# x^{-1} for x > 0.
InvPos <- function(x) { Power(x, -1) }

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

#'
#' The Log class.
#'
#' This class represents the elementwise natural logarithm \eqn{\log(x)}.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name Log-class
#' @aliases Log
#' @rdname Log-class
.Log <- setClass("Log", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname Log-class
Log <- function(x) { .Log(x = x) }

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{Log} object.
#' @param values A list of arguments to the atom.
#' @describeIn Log The elementwise natural logarithm of the input value.
setMethod("to_numeric", "Log", function(object, values) { log(values[[1]]) })

#' @describeIn Log The sign of the atom is unknown.
setMethod("sign_from_args", "Log", function(object) { c(FALSE, FALSE) })

#' @describeIn Log The atom is not convex.
setMethod("is_atom_convex", "Log", function(object) { FALSE })

#' @describeIn Log The atom is concave.
setMethod("is_atom_concave", "Log", function(object) { TRUE })

#' @describeIn Log Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Log", function(object) { FALSE })

#' @describeIn Log Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Log", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn Log The atom is weakly increasing.
setMethod("is_incr", "Log", function(object, idx) { TRUE })

#' @describeIn Log The atom is not weakly decreasing.
setMethod("is_decr", "Log", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Log Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Log", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  # Outside domain or on boundary
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- 1.0/values[[1]]
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

#' @describeIn Log Returns constraints describng the domain of the node
setMethod(".domain", "Log", function(object) { list(object@args[[1]] >= 0) })

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

#'
#' The Logistic class.
#'
#' This class represents the elementwise operation \eqn{\log(1 + e^x)}.
#' This is a special case of log(sum(exp)) that evaluates to a vector rather than to a scalar,
#' which is useful for logistic regression.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name Logistic-class
#' @aliases Logistic
#' @rdname Logistic-class
.Logistic <- setClass("Logistic", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname Logistic-class
Logistic <- function(x) { .Logistic(x = x) }

setMethod("initialize", "Logistic", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{Logistic} object.
#' @param values A list of arguments to the atom.
#' @describeIn Logistic Evaluates \code{e^x} elementwise, adds one, and takes the natural logarithm.
setMethod("to_numeric", "Logistic", function(object, values) { log(1 + exp(values[[1]])) })

#' @describeIn Logistic The atom is positive.
setMethod("sign_from_args", "Logistic", function(object) { c(TRUE, FALSE) })

#' @describeIn Logistic The atom is convex.
setMethod("is_atom_convex", "Logistic", function(object) { TRUE })

#' @describeIn Logistic The atom is not concave.
setMethod("is_atom_concave", "Logistic", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Logistic The atom is weakly increasing.
setMethod("is_incr", "Logistic", function(object, idx) { TRUE })

#' @describeIn Logistic The atom is not weakly decreasing.
setMethod("is_decr", "Logistic", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Logistic Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Logistic", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)
  exp_val <- exp(values[[1]])
  grad_vals <- exp_val/(1 + exp_val)
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

#'
#' The MaxElemwise class.
#'
#' This class represents the elementwise maximum.
#'
#' @slot arg1 The first \linkS4class{Expression} in the maximum operation.
#' @slot arg2 The second \linkS4class{Expression} in the maximum operation.
#' @slot ... Additional \linkS4class{Expression} objects in the maximum operation.
#' @name MaxElemwise-class
#' @aliases MaxElemwise
#' @rdname MaxElemwise-class
.MaxElemwise <- setClass("MaxElemwise", validity = function(object) {
                           if(is.null(object@args) || length(object@args) < 2)
                             stop("[MaxElemwise: validation] args must have at least 2 arguments")
                           return(TRUE)
                         }, contains = "Elementwise")

#' @param arg1 The first \linkS4class{Expression} in the maximum operation.
#' @param arg2 The second \linkS4class{Expression} in the maximum operation.
#' @param ... Additional \linkS4class{Expression} objects in the maximum operation.
#' @rdname MaxElemwise-class
MaxElemwise <- function(arg1, arg2, ...) { .MaxElemwise(atom_args = list(arg1, arg2, ...)) }

#' @param object A \linkS4class{MaxElemwise} object.
#' @param values A list of arguments to the atom.
#' @describeIn MaxElemwise The elementwise maximum.
setMethod("to_numeric", "MaxElemwise", function(object, values) {
  # Reduce(function(x, y) { ifelse(x >= y, x, y) }, values)
  Reduce("pmax", values)
})

#' @describeIn MaxElemwise The sign of the atom.
setMethod("sign_from_args", "MaxElemwise", function(object) {
  # Reduces the list of argument signs according to the following rules:
  #    NONNEGATIVE, ANYTHING = NONNEGATIVE
  #    ZERO, UNKNOWN = NONNEGATIVE
  #    ZERO, ZERO = ZERO
  #    ZERO, NONPOSITIVE = ZERO
  #    UNKNOWN, NONPOSITIVE = UNKNOWN
  #    NONPOSITIVE, NONPOSITIVE = NONPOSITIVE
  is_pos <- any(sapply(object@args, is_nonneg))
  is_neg <- all(sapply(object@args, is_nonpos))
  c(is_pos, is_neg)
})

#' @describeIn MaxElemwise The atom is convex.
setMethod("is_atom_convex", "MaxElemwise", function(object) { TRUE })

#' @describeIn MaxElemwise The atom is not concave.
setMethod("is_atom_concave", "MaxElemwise", function(object) { FALSE })

#' @describeIn MaxElemwise Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MaxElemwise", function(object) { TRUE })

#' @describeIn MaxElemwise Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MaxElemwise", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MaxElemwise The atom is weakly increasing.
setMethod("is_incr", "MaxElemwise", function(object, idx) { TRUE })

#' @describeIn MaxElemwise The atom is not weakly decreasing.
setMethod("is_decr", "MaxElemwise", function(object, idx) { FALSE })

#' @describeIn MaxElemwise Are all the arguments piecewise linear?
setMethod("is_pwl", "MaxElemwise", function(object) { all(sapply(object@args, is_pwl)) })

#' @param values A list of numeric values for the arguments
#' @describeIn MaxElemwise Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MaxElemwise", function(object, values) {
  max_vals <- to_numeric(object, values)
  vals_dim <- dim(max_vals)
  if(is.null(vals_dim))
    unused <- matrix(TRUE, nrow = length(max_vals), ncol = 1)
  else
    unused <- array(TRUE, dim = vals_dim)
  grad_list <- list()
  idx <- 1
  for(value in values) {
    rows <- size(object@args[[idx]])
    cols <- size(object)
    grad_vals <- (value == max_vals) & unused

    # Remove all the max_vals that were used
    unused[value == max_vals] <- FALSE
    grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
    idx <- idx + 1
  }
  grad_list
})

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
#' 
#' An alias for -MinElemwise(x, 0)
#' 
#' @return An alias for -MinElemwise(x, 0)
Neg <- function(x) { -MinElemwise(x, 0) }
#' 
#' An alias for MaxElemwise(x, 0)
#' 
#' @return An alias for MaxElemwise(x, 0)
Pos <- function(x) { MaxElemwise(x, 0) }

#'
#' The Power class.
#'
#' This class represents the elementwise power function \eqn{f(x) = x^p}.
#' If \code{expr} is a CVXR expression, then \code{expr^p} is equivalent to \code{Power(expr, p)}.
#'
#' #' For \eqn{p = 0}, \eqn{f(x) = 1}, constant, positive.
#' For \eqn{p = 1}, \eqn{f(x) = x}, affine, increasing, same sign as \eqn{x}.
#' For \eqn{p = 2,4,8,...}, \eqn{f(x) = |x|^p}, convex, signed monotonicity, positive.
#' For \eqn{p < 0} and \eqn{f(x) = }
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x > 0}}
#'   \item{\eqn{+\infty}}{\eqn{x \leq 0}}
#' }, this function is convex, decreasing, and positive.
#' For \eqn{0 < p < 1} and \eqn{f(x) =}
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x \geq 0}}
#'   \item{\eqn{-\infty}}{\eqn{x < 0}}
#' }, this function is concave, increasing, and positive.
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

#' @describeIn Power Returns constraints describng the domain of the node
setMethod(".domain", "Power", function(object) {
  if((object@p < 1 && object@p != 0) || (object@p > 1 && !is_power2(object@p)))
    list(object@args[[1]] >= 0)
  else
    list()
})

#' @describeIn Power A list containing the output of \code{pow_low, pow_mid}, or \code{pow_high} depending on the input power.
setMethod("get_data", "Power", function(object) { list(object@p, object@w) })

#' @param args A list of arguments to reconstruct the atom. If args=NULL, use the current args of the atom
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

setMethod("name", "Power", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ", ", x@p, ")", sep = "")
})

Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }
