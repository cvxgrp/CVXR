#'
#' The Elementwise class.
#'
#' This virtual class represents an elementwise atom.
#'
#' @aliases Elementwise
#' @export
Elementwise <- setClass("Elementwise", contains = c("VIRTUAL", "Atom"))

setMethod("validate_args", "Elementwise", function(object) {
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
})

setMethod("size_from_args", "Elementwise", function(object) {
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
})

Elementwise.elemwise_grad_to_diag <- function(value, rows, cols) {
  value <- as.numeric(value)
  sparseMatrix(i = 1:rows, j = 1:cols, x = value, dims = c(rows, cols))
}

Elementwise.promote <- function(arg, size) {
  if(any(size(arg) != size))
    lo.promote(arg, size)
  else
    arg
}

#'
#' The Abs class.
#'
#' This class represents the elementwise absolute value.
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Abs
#' @export
.Abs <- setClass("Abs", representation(x = "Expression"), contains = "Elementwise")
Abs <- function(x) { .Abs(x = x) }

setMethod("initialize", "Abs", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("to_numeric", "Abs", function(object, values) { abs(values[[1]]) })
setMethod("sign_from_args", "Abs", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Abs", function(object) { TRUE })
setMethod("is_atom_concave", "Abs", function(object) { FALSE })
setMethod("is_incr", "Abs", function(object, idx) { is_positive(object@args[[idx]]) })
setMethod("is_decr", "Abs", function(object, idx) { is_negative(object@args[[idx]]) })
setMethod("is_pwl", "Abs", function(object) { is_pwl(object@args[[1]]) })

setMethod(".grad", "Abs", function(object, values) {
  # Grad: +1 if positive, -1 if negative
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  
  arg_size <- size(object@args[[1]])
  D <- matrix(0, nrow = arg_size[1], ncol = arg_size[2])
  D <- D + (values[[1]] > 0)
  D <- D - (values[[1]] < 0)
  list(Elementwise.elemwise_grad_to_diag(D, rows, cols))
})

Abs.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size(x))
  constraints <- list(create_geq(lo.sum_expr(list(x, t))), create_leq(x, t))
  list(t, constraints)
}

setMethod("graph_implementation", "Abs", function(object, arg_objs, size, data = NA_real_) {
  Abs.graph_implementation(arg_objs, size, data)
})

#'
#' The Entr class.
#'
#' This class represents the elementwise operation -x * log(x).
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Entr
#' @export
.Entr <- setClass("Entr", representation(x = "ConstValORExpr"), contains = "Elementwise")
Entr <- function(x) { .Entr(x = x) }

setMethod("initialize", "Entr", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

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
  results
})

setMethod("sign_from_args", "Entr", function(object) { c(FALSE, FALSE) })
setMethod("is_atom_convex", "Entr", function(object) { FALSE })
setMethod("is_atom_concave", "Entr", function(object) { TRUE })
setMethod("is_incr", "Entr", function(object, idx) { FALSE })
setMethod("is_decr", "Entr", function(object, idx) { FALSE })

setMethod(".grad", "Entr", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  
  # Outside domain or on boundary
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- -log(values[[1]]) - 1
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

setMethod(".domain", "Entr", function(object) { list(object@args[[1]] >= 0) })

Entr.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  list(t, list(ExpCone(t, x, ones)))
}

setMethod("graph_implementation", "Entr", function(object, arg_objs, size, data = NA_real_) {
  Entr.graph_implementation(arg_objs, size, data)
})

#'
#' The Exp class.
#'
#' This class represents the elementwise natural exponential e^x.
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Exp
#' @export
.Exp <- setClass("Exp", representation(x = "Expression"), contains = "Elementwise")
Exp <- function(x) { .Exp(x = x) }

setMethod("initialize", "Exp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("to_numeric", "Exp", function(object, values) { exp(values[[1]]) })
setMethod("sign_from_args", "Exp", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Exp", function(object) { TRUE })
setMethod("is_atom_concave", "Exp", function(object) { FALSE })
setMethod("is_incr", "Exp", function(object, idx) { TRUE })
setMethod("is_decr", "Exp", function(object, idx) { FALSE })

setMethod(".grad", "Exp", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  grad_vals <- exp(values[[1]])
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

Exp.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  list(t, list(ExpCone(x, ones, t)))
}

setMethod("graph_implementation", "Exp", function(object, arg_objs, size, data = NA_real_) {
  Exp.graph_implementation(arg_objs, size, data)
})

#'
#' The Huber class.
#'
#' This class represents the elementwise Huber function.
#' Huber(x, M) = 2M|x|-M^2 for |x| >= |M| and |x|^2 for |x| <= M
#' M defaults to 1.
#'
#' @slot x A \S4class{Expression} that represents the first argument.
#' @slot M A numeric value that represents the second argument.
#' @aliases Huber
#' @export
.Huber <- setClass("Huber", representation(x = "ConstValORExpr", M = "ConstValORExpr"), 
                           prototype(M = 1), contains = "Elementwise")
Huber <- function(x, M = 1) { .Huber(x = x, M = M) }

setMethod("initialize", "Huber", function(.Object, ..., x, M = 1) {
  .Object@M <- as.Constant(M)
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("validate_args", "Huber", function(object) {
  if(!(is_positive(object@M) && is_constant(object@M) && is_scalar(object@M)))
    stop("M must be a non-negative scalar constant")
})

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
    2*huber_loss(M_val, val)
  else if(is.vector(val))
    2*sapply(val, huber_loss(M_val, v))
  else
    2*apply(val, c(1,2), function(v) { huber_loss(M_val, v) })
})

setMethod("sign_from_args", "Huber", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Huber", function(object) { TRUE })
setMethod("is_atom_concave", "Huber", function(object) { FALSE })
setMethod("is_incr", "Huber", function(object, idx) { is_positive(object@args[[idx]]) })
setMethod("is_decr", "Huber", function(object, idx) { is_negative(object@args[[idx]]) })
setMethod("get_data", "Huber", function(object) { list(object@M) })

setMethod(".grad", "Huber", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))

  val_abs <- abs(values[[1]])
  M_val <- as.numeric(value(object@M))
  min_val <- ifelse(val_abs >= M_val, M_val, val_abs)
  
  grad_vals <- 2*(sign(values[[1]]) * min_val)
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

Huber.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  M <- data[[1]]
  x <- arg_objs[[1]]
  n <- create_var(size)
  s <- create_var(size)
  two <- create_const(2, c(1,1))
  
  if(is(M, "Parameter"))
    M <- create_param(M, c(1,1))
  else   # M is constant
    M <- create_const(value(M), c(1,1))
  
  # n^2 + 2*M*|s|
  power_graph <- Power.graph_implementation(list(n), size, list(2, c(as.bigq(1,2), as.bigq(1,2)) ))  # TODO: Check last argument is correct nested list
  n2 <- power_graph[[1]]
  constr_sq <- power_graph[[2]]
  abs_graph <- Abs.graph_implementation(list(s), size)
  abs_s <- abs_graph[[1]]
  constr_abs <- abs_graph[[2]]
  M_abs_s <- lo.mul_expr(M, abs_s, size)
  obj <- lo.sum_expr(list(n2, lo.mul_expr(two, M_abs_s, size)))
  
  # x == s + n
  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, list(create_eq(x, lo.sum_expr(list(n, s)))))
  list(obj, constraints)
}

setMethod("graph_implementation", "Huber", function(object, arg_objs, size, data = NA_real_) {
  Huber.graph_implementation(arg_objs, size, data)
})

InvPos <- function(x) { Power(x, -1) }

.KLDiv <- setClass("KLDiv", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Elementwise")
KLDiv <- function(x, y) { .KLDiv(x = x, y = y) }

setMethod("initialize", "KLDiv", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y))
})

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

setMethod("sign_from_args", "KLDiv", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "KLDiv", function(object) { TRUE })
setMethod("is_atom_concave", "KLDiv", function(object) { FALSE })
setMethod("is_incr", "KLDiv", function(object, idx) { FALSE })
setMethod("is_decr", "KLDiv", function(object, idx) { FALSE })

setMethod(".grad", "KLDiv", function(object, values) {
  if(min(values[[1]]) <= 0 || min(values[[2]]) <= 0)
    return(list(NA_real_, NA_real_))   # Non-differentiable
  else {
    div <- values[[1]]/values[[2]]
    grad_vals <- list(log(div), 1-div)
    grad_list <- list()
    for(idx in 1:length(values)) {
      rows <- prod(size(object@args[[idx]]))
      cols <- prod(size(object))
      grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals[[idx]], rows, cols)))
    }
    return(grad_list)
  }
})

setMethod(".domain", "KLDiv", function(object) { list(object@args[[1]] >= 0, object@args[[2]] >= 0) })

KLDiv.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- Elementwise.promote(arg_objs[[1]], size)
  y <- Elementwise.promote(arg_objs[[2]], size)
  t <- create_var(size)
  constraints <- list(ExpCone(t, x, y), create_geq(y))   # y >= 0
  
  # -t - x + y
  obj <- lo.sub_expr(y, lo.sum_expr(list(x, t)))
  list(obj, constraints)
}

setMethod("graph_implementation", "KLDiv", function(object, arg_objs, size, data = NA_real_) {
  KLDiv.graph_implementation(arg_objs, size, data)
})

#'
#' The Log class.
#'
#' This class represents the elementwise natural logarithm log(x).
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Log
#' @export
.Log <- setClass("Log", representation(x = "ConstValORExpr"), contains = "Elementwise")
Log <- function(x) { .Log(x = x) }

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("to_numeric", "Log", function(object, values) { log(values[[1]]) })
setMethod("sign_from_args", "Log", function(object) { c(FALSE, FALSE) })
setMethod("is_atom_convex", "Log", function(object) { FALSE })
setMethod("is_atom_concave", "Log", function(object) { TRUE })
setMethod("is_incr", "Log", function(object, idx) { TRUE })
setMethod("is_decr", "Log", function(object, idx) { FALSE })

setMethod(".grad", "Log", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  
  # Outside domain or on boundary
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- 1.0/values[[1]]
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

setMethod(".domain", "Log", function(object) { list(object@args[[1]] >= 0) })

Log.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  list(t, list(ExpCone(t, ones, x)))
}

setMethod("graph_implementation", "Log", function(object, arg_objs, size, data = NA_real_) {
  Log.graph_implementation(arg_objs, size, data)
})

#'
#' The Log1p class.
#'
#' This class represents the elementwise operation log(1 + x).
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Log1p
#' @export
.Log1p <- setClass("Log1p", contains = "Log")
Log1p <- function(x) { .Log1p(x = x) }

setMethod("to_numeric", "Log1p", function(object, values) { log(1+values[[1]]) })
setMethod("sign_from_args", "Log1p", function(object) { c(is_positive(object@args[[1]]), is_negative(object@args[[1]])) })

setMethod(".grad", "Log1p", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  
  # Outside domain or on boundary
  if(min(values[[1]]) <= -1)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- 1.0/(values[[1]] + 1)
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

setMethod(".domain", "Log1p", function(object) { list(object@args[[1]] >= -1) })

Log1p.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  ones <- create_const(matrix(1, nrow = x$size[1], ncol = x$size[2]), x$size)
  xp1 <- lo.sum_expr(list(x, ones))
  Log.graph_implementation(list(xp1), size, data)
}

setMethod("graph_implementation", "Log1p", function(object, arg_objs, size, data = NA_real_) {
  Log1p.graph_implementation(arg_objs, size, data)
})

#'
#' The Logistic class.
#'
#' This class represents the elementwise operation log(1 + e^x).
#' This is a special case of log(sum(exp)) that evaluates to a vector rather than to a scalar,
#' which is useful for logistic regression.
#'
#' @slot x The \S4class{Expression} that is being operated upon.
#' @aliases Logistic
#' @export
.Logistic <- setClass("Logistic", representation(x = "Expression"), contains = "Elementwise")
Logistic <- function(x) { .Logistic(x = x) }

setMethod("initialize", "Logistic", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("to_numeric", "Logistic", function(object, values) { log(1 + exp(values[[1]])) })
setMethod("sign_from_args", "Logistic", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Logistic", function(object) { TRUE })
setMethod("is_atom_concave", "Logistic", function(object) { FALSE })
setMethod("is_incr", "Logistic", function(object, idx) { TRUE })
setMethod("is_decr", "Logistic", function(object, idx) { FALSE })

setMethod(".grad", "Logistic", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  exp_val <- exp(values[[1]])
  grad_vals <- exp_val/(1 + exp_val)
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

Logistic.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size)
  
  # log(1 + exp(x)) <= t <=> exp(-t) + exp(x - t) <= 1
  graph0 <- Exp.graph_implementation(list(lo.neg_expr(t)), size)
  obj0 <- graph0[[1]]
  constr0 <- graph0[[2]]
  
  graph1 <- Exp.graph_implementation(list(lo.sub_expr(x, t)), size)
  obj1 <- graph1[[1]]
  constr1 <- graph1[[2]]
  
  lhs <- lo.sum_expr(list(obj0, obj1))
  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  constr <- c(constr0, constr1, list(create_leq(lhs, ones)))
  list(t, constr)
}

setMethod("graph_implementation", "Logistic", function(object, arg_objs, size, data = NA_real_) {
  Logistic.graph_implementation(arg_objs, size, data)
})

#'
#' The MaxElemwise class.
#'
#' This class represents the elementwise maximum.
#'
#' @slot arg1 The first \S4class{Expression} in the maximum operation.
#' @slot arg2 The second \S4class{Expression} in the maximum operation.
#' @slot ... Additional \S4class{Expression}s in the maximum operation.
#' @aliases MaxElemwise
#' @export
.MaxElemwise <- setClass("MaxElemwise", validity = function(object) {
                           if(is.null(object@args) || length(object@args) < 2)
                             stop("[MaxElemwise: validation] args must have at least 2 arguments")
                           return(TRUE)
                         }, contains = "Elementwise")
MaxElemwise <- function(arg1, arg2, ...) { .MaxElemwise(args = list(arg1, arg2, ...)) }

setMethod("to_numeric", "MaxElemwise", function(object, values) {
  # Reduce(function(x, y) { ifelse(x >= y, x, y) }, values)
  Reduce("pmax", values)
})

setMethod("sign_from_args", "MaxElemwise", function(object) {
  is_pos <- any(sapply(object@args, function(arg) { is_positive(arg) }))
  is_neg <- all(sapply(object@args, function(arg) { is_negative(arg) }))
  c(is_pos, is_neg)
})

setMethod("is_atom_convex", "MaxElemwise", function(object) { TRUE })
setMethod("is_atom_concave", "MaxElemwise", function(object) { FALSE })
setMethod("is_incr", "MaxElemwise", function(object, idx) { TRUE })
setMethod("is_decr", "MaxElemwise", function(object, idx) { FALSE })
setMethod("is_pwl", "MaxElemwise", function(object) { all(sapply(object@args, function(arg) { is_pwl(arg) })) })

setMethod(".grad", "MaxElemwise", function(object, values) {
  max_vals <- to_numeric(object, values)
  dims <- dim(max_vals)
  if(is.null(dims))
    unused <- matrix(TRUE, nrow = length(max_vals), ncol = 1)
  else
    unused <- matrix(TRUE, nrow = dims[1], ncol = dims[2])
  grad_list <- list()
  idx <- 1
  for(value in values) {
    rows <- prod(size(object@args[[idx]]))
    cols <- prod(size(object))
    grad_vals <- (value == max_vals) & unused
    
    # Remove all the max_vals that were used
    unused[value == max_vals] <- FALSE
    grad_list <- c(grad_list, list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
  grad_list
})

MaxElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  constraints <- lapply(arg_objs, function(obj) {
    obj <- Elementwise.promote(obj, size)
    create_leq(obj, t) 
  })
  list(t, constraints)
}

setMethod("graph_implementation", "MaxElemwise", function(object, arg_objs, size, data = NA_real_) {
  MaxElemwise.graph_implementation(arg_objs, size, data)
})

MinElemwise <- function(arg1, arg2, ...) {
  min_args <- lapply(c(list(arg1), list(arg2), list(...)), function(arg) { -as.Constant(arg) })
  -.MaxElemwise(args = min_args)
}

Neg <- function(x) { -MinElemwise(x, 0) }

Pos <- function(x) { MaxElemwise(x, 0) }

#'
#' The Power class.
#'
#' This class represents the elementwise power function f(x) = x^p.
#'
#' @slot x The \S4class{Expression} to be raised to a power.
#' @slot p A numeric value indicating the scalar power.
#' @slot max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @aliases Power
#' @export
.Power <- setClass("Power", representation(x = "ConstValORExpr", p = "NumORgmp", max_denom = "numeric", w = "NumORgmp", approx_error = "numeric"), 
                          prototype(max_denom = 1024, w = NA_real_, approx_error = NA_real_), contains = "Elementwise")

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
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("validate_args", "Power", function(object) { })
setMethod("get_data", "Power", function(object) { list(object@p, object@w) })
setMethod("to_numeric", "Power", function(object, values) {
  # Throw error if negative and Power doesn't handle that
  if(object@p < 0 && min(values[[1]]) <= 0)
    stop("Power cannot be applied to negative or zero values")
  else if(is_power2(object@p) && object@p != 0 && min(values[[1]]) < 0)
    stop("Power cannot be applied to negative values")
  else
    return(values[[1]]^(as.double(object@p)))
})

setMethod("sign_from_args", "Power", function(object) {
  if(object@p == 1)   # Same as input
    c(is_positive(object@args[[1]]), is_negative(object@args[[1]]))
  else   # Always positive
    c(TRUE, FALSE)
})

setMethod("is_atom_convex", "Power", function(object) { object@p <= 0 || object@p >= 1 })
setMethod("is_atom_concave", "Power", function(object) { object@p >= 0 && object@p <= 1 })
setMethod("is_constant", "Power", function(object) { object@p == 0 || callNextMethod() })

setMethod("is_incr", "Power", function(object, idx) {
  if(object@p >= 0 && object@p <= 1)
    return(TRUE)
  else if(object@p > 1) {
    if(is_power2(object@p))
      return(is_positive(object@args[[idx]]))
    else
      return(TRUE)
  } else
    return(FALSE)
})

setMethod("is_decr", "Power", function(object, idx) {
  if(object@p <= 0)
    return(TRUE)
  else if(object@p > 1) {
    if(is_power2(object@p))
      return(is_negative(object@args[[idx]]))
    else
      return(FALSE)
  } else
    return(FALSE)
})

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

setMethod(".grad", "Power", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  
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

setMethod(".domain", "Power", function(object) {
  if((object@p < 1 && object@p != 0) || (object@p > 1 && !is_power2(object@p)))
    list(object@args[[1]] >= 0)
  else
    list()
})

Power.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  p <- data[[1]]
  w <- data[[2]]
  
  if(p == 1)
    return(list(x, list()))
  else {
    one <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
    if(p == 0)
      return(one, list())
    else {
      t <- create_var(size)
      
      # TODO: Temporary hack for powers of 1/2 and 2 until gm_constrs works
      # if(p == 1/2)
      #   return(list(t, gm_constrs_spec(t, list(x, one), w)))
      # else if(p == 2)
      #   return(list(t, gm_constrs_spec(t, list(x, one), w)))
      
      if(p > 0 && p < 1)
        return(list(t, gm_constrs(t, list(x, one), w)))
      else if(p > 1)
        return(list(t, gm_constrs(x, list(t, one), w)))
      else if(p < 0)
        return(list(t, gm_constrs(one, list(x, t), w)))
      else
        stop("This power is not yet supported")
    }
  }
}

setMethod("graph_implementation", "Power", function(object, arg_objs, size, data = NA_real_) {
  Power.graph_implementation(arg_objs, size, data)
})

Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }

# Sqrt <- function(x) { Power(x, 1/2) }
# TODO: Get rid of Sqrt class once Fraction handling is implemented in Power
.Sqrt <- setClass("Sqrt", contains = "Elementwise")
Sqrt <- function(x) { .Sqrt(args = list(x)) }

setMethod("validate_args", "Sqrt", function(object) {})
setMethod("to_numeric", "Sqrt", function(object, values) { values[[1]]^0.5 })
setMethod("get_data", "Sqrt", function(object) { list(0.5, c(0.5, 0.5)) })
setMethod("sign_from_args", "Sqrt", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Sqrt", function(object) { FALSE })
setMethod("is_atom_concave", "Sqrt", function(object) { TRUE })
setMethod("is_incr", "Sqrt", function(object, idx) { TRUE })
setMethod("is_decr", "Sqrt", function(object, idx) { FALSE })
setMethod("is_quadratic", "Sqrt", function(object) { is_constant(object@args[[1]]) })

setMethod(".grad", "Sqrt", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))
  grad_vals <- 0.5*values[[1]]^(-0.5)
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

setMethod(".domain", "Sqrt", function(object) { list(object@args[[1]] >= 0) })

Sqrt.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size)
  one <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  two <- create_const(2, c(1, 1))
  length <- prod(size(x))
  constraints <- list(SOCAxis(lo.reshape(lo.sum_expr(list(x, one)), c(length, 1)),
                              lo.vstack(list(
                              lo.reshape(lo.sub_expr(x, one), c(1, length)),
                              lo.reshape(lo.mul_expr(two, t, size(t)), c(1, length))
                              ), c(2, length)),
                            2))
  list(t, constraints)
}

setMethod("graph_implementation", "Sqrt", function(object, arg_objs, size, data = NA_real_) {
  Sqrt.graph_implementation(arg_objs, size, data)
})

# Square <- function(x) { Power(x, 2) }
# TODO: Get rid of Square class once Fraction object is implemented in Power
.Square <- setClass("Square", contains = "Elementwise")
Square <- function(x) { .Square(args = list(x)) }

setMethod("validate_args", "Square", function(object) {})
setMethod("to_numeric", "Square", function(object, values) { values[[1]]^2 })
setMethod("get_data", "Square", function(object) { list(0.5, c(2,-1)) })
setMethod("sign_from_args", "Square", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Square", function(object) { TRUE })
setMethod("is_atom_concave", "Square", function(object) { FALSE })
setMethod("is_incr", "Square", function(object, idx) { is_positive(object@args[[idx]]) })
setMethod("is_decr", "Square", function(object, idx) { is_negative(object@args[[idx]]) })
setMethod("is_quadratic", "Square", function(object) { is_affine(object@args[[1]]) })

setMethod(".grad", "Square", function(object, values) {
  rows <- prod(size(object@args[[1]]))
  cols <- prod(size(object))
  grad_vals <- 2*values[[1]]
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})

setMethod(".domain", "Square", function(object) { list() })

Square.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size)
  one <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  two <- create_const(2, c(1, 1))
  length <- prod(size(x))
  constraints <- list(SOCAxis(lo.reshape(lo.sum_expr(list(t, one)), c(length, 1)),
                              lo.vstack(list(
                                lo.reshape(lo.sub_expr(t, one), c(1, length)),
                                lo.reshape(lo.mul_expr(two, x, size(x)), c(1, length))
                                ), c(2, length)),
                              2))
  list(t, constraints)
}

setMethod("graph_implementation", "Square", function(object, arg_objs, size, data = NA_real_) {
  Square.graph_implementation(arg_objs, size, data)
})

