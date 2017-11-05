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
#' @describeIn Elementwise Check all the shapes are the same or can be promoted.
setMethod("validate_args", "Elementwise", function(object) {
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
})

#' @describeIn Elementwise Size is the same as the sum of the arguments' sizes.
setMethod("size_from_args", "Elementwise", function(object) {
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
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
  value <- as.numeric(value)
  sparseMatrix(i = 1:rows, j = 1:cols, x = value, dims = c(rows, cols))
}

#
# Promotes LinOp
#
# Promotes the LinOp if necessary.
# @param arg The LinOp to promote.
# @param size The desired size.
# @return The promoted LinOp.
# @rdname Elementwise-promote
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
#' @slot x An \linkS4class{Expression} object.
#' @name Abs-class
#' @aliases Abs
#' @rdname Abs-class
.Abs <- setClass("Abs", representation(x = "Expression"), contains = "Elementwise")

#' @param x,object An \linkS4class{Expression} object.
#' @rdname Abs-class
Abs <- function(x) { .Abs(x = x) }

setMethod("initialize", "Abs", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @param values A list of arguments to the atom.
#' @describeIn Abs The elementwise absolute value of the input.
setMethod("to_numeric", "Abs", function(object, values) { abs(values[[1]]) })

#' @describeIn Abs The atom is positive.
setMethod("sign_from_args", "Abs", function(object) { c(TRUE, FALSE) })

#' @describeIn Abs The atom is convex.
setMethod("is_atom_convex", "Abs", function(object) { TRUE })

#' @describeIn Abs The atom is not concave.
setMethod("is_atom_concave", "Abs", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Abs A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Abs", function(object, idx) { is_positive(object@args[[idx]]) })

#' @describeIn Abs A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Abs", function(object, idx) { is_negative(object@args[[idx]]) })

#' @describeIn Abs Is \code{x} piecewise linear?
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

#' @describeIn Abs The graph implementation of the atom.
setMethod("graph_implementation", "Abs", function(object, arg_objs, size, data = NA_real_) {
  Abs.graph_implementation(arg_objs, size, data)
})

#'
#' The Entr class.
#'
#' This class represents the elementwise operation \eqn{-xlog(x)}.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Entr-class
#' @aliases Entr
#' @rdname Entr-class
.Entr <- setClass("Entr", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @rdname Entr-class
Entr <- function(x) { .Entr(x = x) }

setMethod("initialize", "Entr", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

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
  results
})

#' @describeIn Entr The sign of the atom is unknown.
setMethod("sign_from_args", "Entr", function(object) { c(FALSE, FALSE) })

#' @describeIn Entr The atom is not convex.
setMethod("is_atom_convex", "Entr", function(object) { FALSE })

#' @describeIn Entr The atom is concave.
setMethod("is_atom_concave", "Entr", function(object) { TRUE })

#' @describeIn Entr The atom is weakly increasing.
setMethod("is_incr", "Entr", function(object, idx) { FALSE })

#' @describeIn Entr The atom is weakly decreasing.
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

#' @describeIn Entr The graph implementation of the atom.
setMethod("graph_implementation", "Entr", function(object, arg_objs, size, data = NA_real_) {
  Entr.graph_implementation(arg_objs, size, data)
})

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

#' @rdname Exp-class
Exp <- function(x) { .Exp(x = x) }

setMethod("initialize", "Exp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn Exp The matrix with each element exponentiated.
setMethod("to_numeric", "Exp", function(object, values) { exp(values[[1]]) })

#' @describeIn Exp The atom is positive.
setMethod("sign_from_args", "Exp", function(object) { c(TRUE, FALSE) })

#' @describeIn Exp The atom is convex.
setMethod("is_atom_convex", "Exp", function(object) { TRUE })

#' @describeIn Exp The atom is not concave.
setMethod("is_atom_concave", "Exp", function(object) { FALSE })

#' @describeIn Exp The atom is weakly increasing.
setMethod("is_incr", "Exp", function(object, idx) { TRUE })

#' @describeIn Exp The atom is not weakly decreasing.
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

#' @describeIn Exp The graph implementation of the atom.
setMethod("graph_implementation", "Exp", function(object, arg_objs, size, data = NA_real_) {
  Exp.graph_implementation(arg_objs, size, data)
})

#'
#' The Huber class.
#'
#' This class represents the elementwise Huber function.
#' \deqn{Huber(x, M) = \begin{cases}
#'       2M|x|-M^2 & \mbox{for } |x| \geq |M| \\
#'       |x|^2 & \mbox{for } |x| \leq M
#' \end{cases}}
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @slot M A positive scalar value representing the threshold. Defaults to 1.
#' @name Huber-class
#' @aliases Huber
#' @rdname Huber-class
.Huber <- setClass("Huber", representation(x = "ConstValORExpr", M = "ConstValORExpr"),
                           prototype(M = 1), contains = "Elementwise")

#' @rdname Huber-class
Huber <- function(x, M = 1) { .Huber(x = x, M = M) }

setMethod("initialize", "Huber", function(.Object, ..., x, M = 1) {
  .Object@M <- as.Constant(M)
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn Huber Check that \code{M} is a non-negative constant.
setMethod("validate_args", "Huber", function(object) {
  if(!(is_positive(object@M) && is_constant(object@M) && is_scalar(object@M)))
    stop("M must be a non-negative scalar constant")
})

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
    2*huber_loss(M_val, val)
  else if(is.vector(val))
    2*sapply(val, function(v) { huber_loss(M_val, v) })
  else
    2*apply(val, c(1,2), function(v) { huber_loss(M_val, v) })
})

#' @describeIn Huber The atom is positive.
setMethod("sign_from_args", "Huber", function(object) { c(TRUE, FALSE) })

#' @describeIn Huber The atom is convex.
setMethod("is_atom_convex", "Huber", function(object) { TRUE })

#' @describeIn Huber The atom is not concave.
setMethod("is_atom_concave", "Huber", function(object) { FALSE })

#' @describeIn Huber A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Huber", function(object, idx) { is_positive(object@args[[idx]]) })

#' @describeIn Huber A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Huber", function(object, idx) { is_negative(object@args[[idx]]) })

#' @describeIn Huber A list containing the parameter \code{M}.
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

#' @describeIn Huber The graph implementation of the atom.
setMethod("graph_implementation", "Huber", function(object, arg_objs, size, data = NA_real_) {
  Huber.graph_implementation(arg_objs, size, data)
})

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

#' @rdname KLDiv-class
KLDiv <- function(x, y) { .KLDiv(x = x, y = y) }

setMethod("initialize", "KLDiv", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y))
})

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

#' @describeIn KLDiv The atom is not monotonic in any argument.
setMethod("is_incr", "KLDiv", function(object, idx) { FALSE })

#' @describeIn KLDiv The atom is not monotonic in any argument.
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

#' @describeIn KLDiv The graph implementation of the atom.
setMethod("graph_implementation", "KLDiv", function(object, arg_objs, size, data = NA_real_) {
  KLDiv.graph_implementation(arg_objs, size, data)
})

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

#' @rdname Log-class
Log <- function(x) { .Log(x = x) }

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn Log The elementwise natural logarithm of the input value.
setMethod("to_numeric", "Log", function(object, values) { log(values[[1]]) })

#' @describeIn Log The sign of the atom is unknown.
setMethod("sign_from_args", "Log", function(object) { c(FALSE, FALSE) })

#' @describeIn Log The atom is not convex.
setMethod("is_atom_convex", "Log", function(object) { FALSE })

#' @describeIn Log The atom is concave.
setMethod("is_atom_concave", "Log", function(object) { TRUE })

#' @describeIn Log The atom is weakly increasing.
setMethod("is_incr", "Log", function(object, idx) { TRUE })

#' @describeIn Log The atom is not weakly decreasing.
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

#' @describeIn Log The graph implementation of the atom.
setMethod("graph_implementation", "Log", function(object, arg_objs, size, data = NA_real_) {
  Log.graph_implementation(arg_objs, size, data)
})

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

#' @rdname Log1p-class
Log1p <- function(x) { .Log1p(x = x) }

#' @describeIn Log1p The elementwise natural logarithm of one plus the input value.
setMethod("to_numeric", "Log1p", function(object, values) { log(1+values[[1]]) })

#' @describeIn Log1p The sign of the atom.
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

#' @describeIn Log1p The graph implementation of the atom.
setMethod("graph_implementation", "Log1p", function(object, arg_objs, size, data = NA_real_) {
  Log1p.graph_implementation(arg_objs, size, data)
})

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
.Logistic <- setClass("Logistic", representation(x = "Expression"), contains = "Elementwise")

#' @rdname Logistic-class
Logistic <- function(x) { .Logistic(x = x) }

setMethod("initialize", "Logistic", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn Logistic Evaluates \code{e^x} elementwise, adds one, and takes the natural logarithm.
setMethod("to_numeric", "Logistic", function(object, values) { log(1 + exp(values[[1]])) })

#' @describeIn Logistic The atom is positive.
setMethod("sign_from_args", "Logistic", function(object) { c(TRUE, FALSE) })

#' @describeIn Logistic The atom is convex.
setMethod("is_atom_convex", "Logistic", function(object) { TRUE })

#' @describeIn Logistic The atom is not concave.
setMethod("is_atom_concave", "Logistic", function(object) { FALSE })

#' @describeIn Logistic The atom is weakly increasing.
setMethod("is_incr", "Logistic", function(object, idx) { TRUE })

#' @describeIn Logistic The atom is not weakly decreasing.
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

#' @describeIn Logistic The graph implementation of the atom.
setMethod("graph_implementation", "Logistic", function(object, arg_objs, size, data = NA_real_) {
  Logistic.graph_implementation(arg_objs, size, data)
})

#'
#' The MaxElemwise class.
#'
#' This class represents the elementwise maximum.
#'
#' @slot arg1 The first \linkS4class{Expression} in the maximum operation.
#' @slot arg2 The second \linkS4class{Expression} in the maximum operation.
#' @slot ... Additional \linkS4class{Expression}s in the maximum operation.
#' @name MaxElemwise-class
#' @aliases MaxElemwise
#' @rdname MaxElemwise-class
.MaxElemwise <- setClass("MaxElemwise", validity = function(object) {
                           if(is.null(object@args) || length(object@args) < 2)
                             stop("[MaxElemwise: validation] args must have at least 2 arguments")
                           return(TRUE)
                         }, contains = "Elementwise")

#' @rdname MaxElemwise-class
MaxElemwise <- function(arg1, arg2, ...) { .MaxElemwise(args = list(arg1, arg2, ...)) }

#' @describeIn MaxElemwise The elementwise maximum.
setMethod("to_numeric", "MaxElemwise", function(object, values) {
  # Reduce(function(x, y) { ifelse(x >= y, x, y) }, values)
  Reduce("pmax", values)
})

#' @describeIn MaxElemwise The sign of the atom.
setMethod("sign_from_args", "MaxElemwise", function(object) {
  is_pos <- any(sapply(object@args, function(arg) { is_positive(arg) }))
  is_neg <- all(sapply(object@args, function(arg) { is_negative(arg) }))
  c(is_pos, is_neg)
})

#' @describeIn MaxElemwise The atom is convex.
setMethod("is_atom_convex", "MaxElemwise", function(object) { TRUE })

#' @describeIn MaxElemwise The atom is not concave.
setMethod("is_atom_concave", "MaxElemwise", function(object) { FALSE })

#' @describeIn MaxElemwise The atom is weakly increasing.
setMethod("is_incr", "MaxElemwise", function(object, idx) { TRUE })

#' @describeIn MaxElemwise The atom is not weakly decreasing.
setMethod("is_decr", "MaxElemwise", function(object, idx) { FALSE })

#' @describeIn MaxElemwise Are all the arguments piecewise linear?
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

#' @describeIn MaxElemwise The graph implementation of the atom.
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
#' This class represents the elementwise power function \eqn{f(x) = x^p}.
#' If \code{expr} is a CVXR expression, then \code{expr^p} is equivalent to \code{Power(expr, p)}.
#'
#' \deqn{\begin{array}{ccl}
#' p = 0 & f(x) = 1 & \text{constant, positive} \\
#' p = 1 & f(x) = x & \text{affine, increasing, same sign as $x$} \\
#' p = 2,4,8,\ldots & f(x) = |x|^p  & \text{convex, signed monotonicity, positive} \\
#' p < 0 & f(x) = \begin{cases} x^p & x > 0 \\ +\infty & x \leq 0 \end{cases} & \text{convex, decreasing, positive} \\
#' 0 < p < 1 & f(x) = \begin{cases} x^p & x \geq 0 \\ -\infty & x < 0 \end{cases} & \text{concave, increasing, positive} \\
#' p > 1,\ p \neq 2,4,8,\ldots & f(x) = \begin{cases} x^p & x \geq 0 \\ +\infty & x < 0 \end{cases} & \text{convex, increasing, positive}.
#' \end{array}}
#'
#' @slot x The \linkS4class{Expression} to be raised to a power.
#' @slot p A numeric value indicating the scalar power.
#' @slot max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @name Power-class
#' @aliases Power
#' @rdname Power-class
.Power <- setClass("Power", representation(x = "ConstValORExpr", p = "NumORgmp", max_denom = "numeric", w = "NumORgmp", approx_error = "numeric"),
                          prototype(max_denom = 1024, w = NA_real_, approx_error = NA_real_), contains = "Elementwise")

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
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn Power Verification of arguments happens during initialization.
setMethod("validate_args", "Power", function(object) { return() })

#' @describeIn Power A list containing the output of \code{pow_low, pow_mid}, or \code{pow_high} depending on the input power.
setMethod("get_data", "Power", function(object) { list(object@p, object@w) })

#' @describeIn Power Throw an error if the power is negative and cannot be handled.
setMethod("to_numeric", "Power", function(object, values) {
  # Throw error if negative and Power doesn't handle that
  if(object@p < 0 && min(values[[1]]) <= 0)
    stop("Power cannot be applied to negative or zero values")
  else if(is_power2(object@p) && object@p != 0 && min(values[[1]]) < 0)
    stop("Power cannot be applied to negative values")
  else
    return(values[[1]]^(as.double(object@p)))
})

#' @describeIn Power The sign of the atom.
setMethod("sign_from_args", "Power", function(object) {
  if(object@p == 1)   # Same as input
    c(is_positive(object@args[[1]]), is_negative(object@args[[1]]))
  else   # Always positive
    c(TRUE, FALSE)
})

#' @describeIn Power Is \eqn{p \leq 0} or \eqn{p \geq 1}?
setMethod("is_atom_convex", "Power", function(object) { object@p <= 0 || object@p >= 1 })

#' @describeIn Power Is \eqn{p \geq 0} or \eqn{p \leq 1}?
setMethod("is_atom_concave", "Power", function(object) { object@p >= 0 && object@p <= 1 })

#' @describeIn Power A logical value indicating whether the atom is constant.
setMethod("is_constant", "Power", function(object) { object@p == 0 || callNextMethod() })

#' @describeIn Power A logical value indicating whether the atom is weakly increasing.
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

#' @describeIn Power A logical value indicating whether the atom is weakly decreasing.
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

#' @describeIn Power The graph implementation of the atom.
setMethod("graph_implementation", "Power", function(object, arg_objs, size, data = NA_real_) {
  Power.graph_implementation(arg_objs, size, data)
})

Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }

#'
#' The Sqrt class.
#'
#' This class represents the elementwise square root \eqn{\sqrt{x}}.
#' 
#' @slot x An \linkS4class{Expression} object.
#' @name Sqrt-class
#' @aliases Sqrt
#' @rdname Sqrt-class
.Sqrt <- setClass("Sqrt", contains = "Elementwise")

#' @rdname Sqrt-class
Sqrt <- function(x) { .Sqrt(args = list(x)) }
# Sqrt <- function(x) { Power(x, 1/2) }
# TODO: Get rid of Sqrt class once Fraction handling is implemented in Power

#' @describeIn Sqrt Verification of arguments happens during initialization.
setMethod("validate_args", "Sqrt", function(object) { return() })

#' @describeIn Sqrt The elementwise square root of the input value.
setMethod("to_numeric", "Sqrt", function(object, values) { values[[1]]^0.5 })

#' @describeIn Sqrt A list containing the output of \code{pow_mid}.
setMethod("get_data", "Sqrt", function(object) { list(0.5, c(0.5, 0.5)) })

#' @describeIn Sqrt The atom is positive.
setMethod("sign_from_args", "Sqrt", function(object) { c(TRUE, FALSE) })

#' @describeIn Sqrt The atom is not convex.
setMethod("is_atom_convex", "Sqrt", function(object) { FALSE })

#' @describeIn Sqrt The atom is concave.
setMethod("is_atom_concave", "Sqrt", function(object) { TRUE })

#' @describeIn Sqrt The atom is weakly increasing.
setMethod("is_incr", "Sqrt", function(object, idx) { TRUE })

#' @describeIn Sqrt The atom is not weakly decreasing.
setMethod("is_decr", "Sqrt", function(object, idx) { FALSE })

#' @describeIn Sqrt Is \code{x} constant?
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

#' @describeIn Sqrt The graph implementation of the atom.
setMethod("graph_implementation", "Sqrt", function(object, arg_objs, size, data = NA_real_) {
  Sqrt.graph_implementation(arg_objs, size, data)
})

#'
#' The Square class.
#'
#' This class represents the elementwise square \eqn{x^2}.
#' 
#' @slot x An \linkS4class{Expression}.
#' @name Square-class
#' @aliases Square
#' @rdname Square-class
.Square <- setClass("Square", contains = "Elementwise")

#' @rdname Square-class
Square <- function(x) { .Square(args = list(x)) }
# Square <- function(x) { Power(x, 2) }
# TODO: Get rid of Square class once Fraction object is implemented in Power

#' @describeIn Square Verification of arguments happens during initialization.
setMethod("validate_args", "Square", function(object) { return() })

#' @describeIn Square The elementwise square of the input value.
setMethod("to_numeric", "Square", function(object, values) { values[[1]]^2 })

#' @describeIn Square A list containing the output of \code{pow_high}.
setMethod("get_data", "Square", function(object) { list(0.5, c(2,-1)) })

#' @describeIn Square The atom is positive.
setMethod("sign_from_args", "Square", function(object) { c(TRUE, FALSE) })

#' @describeIn Square The atom is convex.
setMethod("is_atom_convex", "Square", function(object) { TRUE })

#' @describeIn Square The atom is not concave.
setMethod("is_atom_concave", "Square", function(object) { FALSE })

#' @describeIn Square A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "Square", function(object, idx) { is_positive(object@args[[idx]]) })

#' @describeIn Square A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "Square", function(object, idx) { is_negative(object@args[[idx]]) })

#' @describeIn Square Is \code{x} affine?
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

#' @describeIn Square The graph implementation of the atom.
setMethod("graph_implementation", "Square", function(object, arg_objs, size, data = NA_real_) {
  Square.graph_implementation(arg_objs, size, data)
})
