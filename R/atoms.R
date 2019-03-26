#'
#' The Atom class.
#'
#' This virtual class represents atomic expressions in CVXR.
#'
#' @name Atom-class
#' @aliases Atom
#' @rdname Atom-class
Atom <- setClass("Atom", representation(atom_args = "list", .dim = "NumORNULL"), prototype(atom_args = list(), .dim = NULL),
                 validity = function(object) {
                   if(length(object@atom_args) == 0)
                     stop("[Atom: atom_args] no arguments given to ", class(object))
                   return(TRUE)
                 }, contains = c("VIRTUAL", "Expression"))

setMethod("initialize", "Atom", function(.Object, ..., atom_args = list(), .dim = NULL) {
  .Object@args <- lapply(atom_args, as.Constant)
  validate_args(.Object)
  .Object@.dim <- dim_from_args(.Object)
  if(length(.Object@.dim) > 2)
    stop("Atoms must be at most 2D.")
  .Object
})

#' @param x,object An \linkS4class{Atom} object.
setMethod("show", "Atom", function(object) {
  if(is.null(get_data(object)))
    data <- list()
  else
    data <- sapply(get_data(object), as.character)
  arg_names <- sapply(object@args, name)
  cat(class(object), "(", paste(c(arg_names, data), collapse = ", "), ")", sep = "")
})

setMethod("name", "Atom", function(x) {
  if(is.null(get_data(x)))
    data <- list()
  else
    data <- sapply(get_data(x), as.character)
  arg_names <- sapply(x@args, name)
  paste(class(x), "(", paste(c(arg_names, data), collapse = ", "), ")", sep = "")
})

#' @describeIn Atom Raises an error if the arguments are invalid.
setMethod("validate_args", "Atom", function(object) { 
  if(!allow_complex(object) && any(sapply(object@args, is_complex)))
    stop("Arguments to ", class(object), " cannot be complex.")
})

#' @rdname dim_from_args
setMethod("dim_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom The \code{c(row, col)} dimensions of the atom.
setMethod("dim", "Atom", function(x) { x@.dim })

#' @describeIn Atom The number of rows in the atom.
setMethod("nrow", "Atom", function(x) { dim(x)[1] })

#' @describeIn Atom The number of columns in the atom.
setMethod("ncol", "Atom", function(x) { dim(x)[2] })

#' @describeIn Atom Does the atom handle complex numbers?
setMethod("allow_complex", "Atom", function(object) { FALSE })

#' @rdname sign_from_args
setMethod("sign_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom A logical value indicating whether the atom is nonnegative.
setMethod("is_nonneg", "Atom", function(object) { sign_from_args(object)[1] })

#' @describeIn Atom A logical value indicating whether the atom is nonpositive.
setMethod("is_nonpos", "Atom", function(object) { sign_from_args(object)[2] })

#' @describeIn Atom A logical value indicating whether the atom is imaginary.
setMethod("is_imag", "Atom", function(object) { FALSE })

#' @describeIn Atom A logical value indicating whether the atom is complex valued.
setMethod("is_complex", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_convex", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_concave", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_affine", "Atom", function(object) { is_atom_concave(object) && is_atom_convex(object) })

#' @rdname curvature-atom
setMethod("is_atom_log_log_convex", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_log_log_concave", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_log_log_affine", "Atom", function(object) { is_atom_log_log_concave(object) && is_atom_log_log_convex(object) })

#' @rdname curvature-comp
setMethod("is_incr", "Atom", function(object, idx) { stop("Unimplemented") })

#' @rdname curvature-comp
setMethod("is_decr", "Atom", function(object, idx) { stop("Unimplemented") })

#' @describeIn Atom A logical value indicating whether the atom is convex.
setMethod("is_convex", "Atom", function(object) {
  # Applies DCP composition rule
  if(is_constant(object))
    return(TRUE)
  else if(is_atom_convex(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_affine(arg) || (is_convex(arg) && is_incr(object, idx)) || (is_concave(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is concave.
setMethod("is_concave", "Atom", function(object) {
  # Applies DCP composition rule
  if(is_constant(object))
    return(TRUE)
  else if(is_atom_concave(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_affine(arg) || (is_concave(arg) && is_incr(object, idx)) || (is_convex(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is log-log convex.
setMethod("is_log_log_convex", "Atom", function(object) {
  # Verifies DGP composition rule.
  if(is_log_log_constant(object))
    return(TRUE)
  else if(is_atom_log_log_convex(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_log_log_affine(arg) || (is_log_log_convex(arg) && is_incr(object, idx)) || (is_log_log_concave(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is log-log concave.
setMethod("is_log_log_concave", "Atom", function(object) {
  # Verifies DGP composition rule.
  if(is_log_log_constant(object))
    return(TRUE)
  else if(is_atom_log_log_concave(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_log_log_affine(arg) || (is_log_log_concave(arg) && is_incr(object, idx)) || (is_log_log_convex(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom Represent the atom as an affine objective and conic constraints.
setMethod("canonicalize", "Atom", function(object) {
  # Constant atoms are treated as a leaf.
  if(is_constant(object)) {
    # Parameterized expressions are evaluated later.
    if(!is.na(parameters(object)) && length(parameters(object)) > 0) {
      param <- CallbackParam(value(object), dim(object))
      return(canonical_form(param))
    # Non-parameterized expressions are evaluated immediately.
    } else
      return(canonical_form(Constant(value(object))))
  } else {
    arg_objs <- list()
    constraints <- list()
    for(arg in object@args) {
      canon <- canonical_form(arg)
      arg_objs[[length(arg_objs) + 1]] <- canon[[1]]
      constraints <- c(constraints, canon[[2]])
    }
    # Special info required by the graph implementation.
    data <- get_data(object)
    graph <- graph_implementation(object, arg_objs, dim(object), data)
    return(list(graph[[1]], c(constraints, graph[[2]])))
  }
})

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Atom The graph implementation of the atom.
setMethod("graph_implementation", "Atom", function(object, arg_objs, dim, data = NA_real_) { stop("Unimplemented") })

# .value_impl.Atom <- function(object) {
setMethod("value_impl", "Atom", function(object) {
  obj_dim <- dim(object)
  # dims with 0's dropped in presolve.
  if(0 %in% obj_dim)
    result <- matrix(nrow = 0, ncol = 0)
  # Catch the case when the expression is known to be zero through DCP analysis
  else if(is_zero(object))
    result <- matrix(0, nrow = obj_dim[1], ncol = obj_dim[2])
  else {
    arg_values <- list()
    for(arg in object@args) {
      # An argument without a value makes all higher level values NA.
      # But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
      arg_val <- value_impl(arg)
      if(any(is.na(arg_val)) && !is_constant(object))
        return(NA_real_)
      else
        arg_values <- c(arg_values, list(arg_val))
    }
    result <- to_numeric(object, arg_values)
  }
  return(result)
})

#' @describeIn Atom The value of the atom.
setMethod("value", "Atom", function(object) {
  if(any(sapply(parameters(object), function(p) { is.na(value(p)) })))
    return(NA_real_)
  return(value_impl(object))
})

#' @describeIn Atom The (sub/super)-gradient of the atom with respect to each variable.
setMethod("grad", "Atom", function(object) {
  # Short-circuit to all zeros if known to be constant
  if(is_constant(object))
    return(constant_grad(object))

  # Returns NA if variable values are not supplied
  arg_values <- list()
  for(arg in object@args) {
    arg_val <- value(arg)
    if(any(is.na(arg_val)))
      return(error_grad(object))
    else
      arg_values <- c(arg_values, list(arg_val))
  }

  # A list of gradients wrt arguments
  grad_self <- .grad(object, arg_values)

  # The chain rule
  result <- list()
  idx <- 1
  for(arg in object@args) {
    # A dictionary of gradients wrt variables.
    # Partial argument / partial x.
    grad_arg <- grad(arg)
    for(key in names(grad_arg)) {
      # None indicates gradient is not defined.
      if(any(is.na( as.numeric(grad_arg[[key]]) )) || any(is.na( as.numeric(grad_self[[idx]]) )))
        result[[key]] <- NA_real_
      else {
        D <- grad_arg[[key]] %*% grad_self[[idx]]
        # Convert 1x1 matrices to scalars.
        if((is.matrix(D) || is(D, "Matrix")) && all(dim(D) == c(1,1)))
          D <- D[1,1]

        if(key %in% names(result))
          result[[key]] <- result[[key]] + D
        else
          result[[key]] <- D
      }
    }
    idx <- idx + 1
  }
  return(result)
})

setMethod(".grad", "Atom", function(object, values) { stop("Unimplemented") })

#' @describeIn Atom A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "Atom", function(object) {
  cons <- list()
  for(arg in object@args) {
    for(con in domain(arg))
      cons <- c(cons, con)
  }
  c(.domain(object), cons)
})

setMethod(".domain", "Atom", function(object) { list() })

setMethod("atoms", "Atom", function(object) {
  atom_list <- list()
  for(arg in object@args)
    atom_list <- c(atom_list, atoms(arg))
  atom_list <- c(atom_list, list(class(object)))
  return(unique(atom_list))
})

#'
#' The AxisAtom class.
#'
#' This virtual class represents atomic expressions that can be applied along an axis in CVXR.
#'
#' @slot expr A numeric element, data.frame, matrix, vector, or Expression.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name AxisAtom-class
#' @aliases AxisAtom
#' @rdname AxisAtom-class
AxisAtom <- setClass("AxisAtom", representation(expr = "ConstValORExpr", axis = "ANY", keepdims = "logical"), 
                                 prototype(axis = NA_real_, keepdims = FALSE), contains = c("VIRTUAL", "Atom"))

setMethod("initialize", "AxisAtom", function(.Object, ..., expr, axis = NA_real_, keepdims = FALSE) {
  .Object@expr <- expr
  .Object@axis <- axis
  .Object@keepdims <- keepdims
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Atom} object.
#' @describeIn AxisAtom The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "AxisAtom", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(object@keepdims && is.na(object@axis))   # Copy scalar to maintain original dimensions.
    arg_dim <- rep(1, length(arg_dim))
  else if(object@keepdims && !is.na(object@axis)) {   # Collapse dimensions NOT in axis to 1.
    collapse <- setdiff(1:length(arg_dim), object@axis)
    arg_dim[collapse] <- 1
  } else if(!object@keepdims && is.na(object@axis))   # Return a scalar.
    arg_dim <- NULL   # TODO: Should this be NA instead?
  else   # Drop dimensions NOT in axis and collapse atom.
    arg_dim <- arg_dim[object@axis]
  return(arg_dim)
})

#' @describeIn AxisAtom A list containing \code{axis} and \code{keepdims}.
setMethod("get_data", "AxisAtom", function(object) { list(object@axis, object@keepdims) })

#' @describeIn AxisAtom Check that the new dimensions have the same number of entries as the old.
setMethod("validate_args", "AxisAtom", function(object) {
  if(!is.na(object@axis) && any(object@axis > ndim(object@args[[1]]) || object@axis <= 0))
    stop("Invalid argument for axis. Must be an integer between 1 and ", ndim(object@args[[1]]))
  callNextMethod()
})

setMethod(".axis_grad", "AxisAtom", function(object, values) {
  if(is.na(object@axis) || ndim(object@args[[1]]) < 2) {
    value <- matrix(values[[1]], nrow = size(object@args[[1]]), ncol = 1)
    D <- .column_grad(object, value)
    if(is(D, "Matrix") || !any(is.na(D)))
      D <- Matrix(D, sparse = TRUE)
  } else {
    m <- nrow(object@args[[1]])
    n <- ncol(object@args[[1]])
    if(object@axis == 2) {   # Function apply to each column
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, n))
      for(i in 1:n) {
        value <- values[[1]][,i]
        d <- t(.column_grad(object, value))
        if(any(is.na(as.numeric(d))))
          return(list(NA_real_))
        row <- seq((i-1)*n+1, (i-1)*n+m, length.out = m)
        col <- rep(1,m) * i
        D <- D + sparseMatrix(i = row, j = col, x = as.numeric(d), dims = c(m*n, n))
      }
    } else {   # Function apply to each row
      values <- t(values[[1]])
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, m))
      for(i in 1:m) {
        value <- values[,i]
        d <- t(.column_grad(object, value))
        if(any(is.na(as.numeric(d))))
          return(list(NA_real_))
        row <- seq(i, i+(n-1)*m, length.out = n)
        col <- rep(1,n)*i
        D <- D + sparseMatrix(i = row, j = col, x = as.numeric(d), dims = c(m*n, m))
      }
    }
  }
  list(D)
})

setMethod(".column_grad", "AxisAtom", function(object, value) { stop("Unimplemented") })

#'
#' The CumMax class.
#' 
#' This class represents the cumulative maximum of an expression.
#' 
#' @slot x An \linkS4class{Expression}.
#' @name CumMax-class
#' @aliases CumMax
#' @rdname CumMax-class
.CumMax <- setClass("CumMax", contains = "AxisAtom")

#' @param x An \linkS4class{Expression}
#' @param axis A numeric vector indicating the axes along which to apply the function. For a 2D matrix, \code{1} indicates rows, \code{2} indicates columns, and \code{c(1,2)} indicates rows and columns.
#' @rdname CumMax-class
CumMax <- function(x, axis = NA_real_) { .CumMax(expr = x, axis = axis) }

#' @param object A \linkS4class{CumMax} object.
#' @param values A list of arguments to the atom.
#' @rdname CumMax-class The cumulative maximum along the axis.
setMethod("to_numeric", "CumMax", function(object, values) { apply(values[[1]], object@axis, max) })

setMethod(".grad", "CumMax", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "CumMax", function(object, value) {
  # Grad: 1 for a largest index.
  value <- as.numeric(value)
  idx <- (value == max(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  return(D)
})

#' @describeIn CumMax The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "CumMax", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn CumMax Is the atom convex?
setMethod("is_atom_convex", "CumMax", function(object) { TRUE })

#' @describeIn CumMax Is the atom concave?
setMethod("is_atom_concave", "CumMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn CumMax Is the atom weakly increasing in the index?
setMethod("is_incr", "CumMax", function(object, idx) { TRUE })

#' @describeIn CumMax Is the atom weakly decreasing in the index?
setMethod("is_decr", "CumMax", function(object, idx) { FALSE })

#'
#' The EyeMinusInv class.
#'
#' This class represents the unity resolvent of an elementwise positive matrix \eqn{X}, i.e., \eqn{(I - X)^{-1}},
#' and it enforces the constraint that the spectral radius of \eqn{X} is at most 1.
#' This atom is log-log convex.
#' 
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @name EyeMinusInv-class
#' @aliases EyeMinusInv
#' @rdname EyeMinusInv-class
.EyeMinusInv <- setClass("EyeMinusInv", representation(X = "ConstValORExpr"), 
                         validity = function(object) {
                           if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                             stop("[EyeMinusInv: X] The argument X must be a square matrix.")
                           return(TRUE)
                          }, contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname EyeMinusInv-class
EyeMinusInv <- function(X) { .EyeMinusInv(X = X) }

setMethod("initialize", "EyeMinusInv", function(.Object, ..., X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param object A \linkS4class{EyeMinusInv} object.
#' @param values A list of arguments to the atom.
#' @describeIn EyeMinusInv The unity resolvent of the matrix.
setMethod("to_numeric", "EyeMinusInv", function(object, values) {
  base::solve(diag(nrow(object@args[[1]])) - values[[1]])
})

#' @param x,object An \linkS4class{EyeMinusInv} object.
setMethod("name", "EyeMinusInv", function(x) { paste(class(x), x@args[[1]]) })

#' @describeIn EyeMinusInv The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "EyeMinusInv", function(object) { dim(object@args[[1]]) })

#' @describeIn EyeMinusInv The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "EyeMinusInv", function(object) { c(TRUE, FALSE) })

#' @describeIn EyeMinusInv Is the atom convex?
setMethod("is_atom_convex", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom concave?
setMethod("is_atom_concave", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "EyeMinusInv", function(object) { TRUE })

#' @describeIn EyeMinusInv Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "EyeMinusInv", function(object) { FALSE })

# TODO: Figure out monotonicity.
#' @param idx An index into the atom.
#' @describeIn EyeMinusInv Is the atom weakly increasing in the index?
setMethod("is_incr", "EyeMinusInv", function(object, idx) { FALSE })

#' @describeIn EyeMinusInv Is the atom weakly decreasing in the index?
setMethod("is_decr", "EyeMinusInv", function(object, idx) { FALSE })

setMethod(".grad", "EyeMinusInv", function(object, values) { NA_real_ })

# The resolvent of a positive matrix, (sI - X)^(-1).
# For an elementwise positive matrix X and a positive scalar s, this atom computes
# (sI - X)^(-1), and it enforces the constraint that the spectral radius of X/s is
# at most 1.
# This atom is log-log convex.
Resolvent <- function(X, s) {
  1.0 / (s * EyeMinusInv(X / s))
}

#'
#' The GeoMean class.
#'
#' This class represents the (weighted) geometric mean of vector \eqn{x} with optional powers given by \eqn{p}.
#'
#' \deqn{\left(x_1^{p_1} \cdots x_n^{p_n} \right)^{\frac{1}{\mathbf{1}^Tp}}}
#'
#' The geometric mean includes an implicit constraint that \eqn{x_i \geq 0} whenever \eqn{p_i > 0}. If \eqn{p_i = 0, x_i} will be unconstrained.
#' The only exception to this rule occurs when \eqn{p} has exactly one nonzero element, say \eqn{p_i}, in which case \code{GeoMean(x,p)} is equivalent to \eqn{x_i} (without the nonnegativity constraint).
#' A specific case of this is when \eqn{x \in \mathbf{R}^1}.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @slot max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @slot w (Internal) A list of \code{bigq} objects that represent a rational approximation of \code{p/sum(p)}.
#' @slot approx_error (Internal) The error in approximating \code{p/sum(p)} with \code{w}, given by \eqn{\|p/\mathbf{1}^Tp - w\|_{\infty}}.
#' @name GeoMean-class
#' @aliases GeoMean
#' @importClassesFrom gmp bigq bigz
#' @rdname GeoMean-class
.GeoMean <- setClass("GeoMean", representation(x = "ConstValORExpr", p = "numeric", max_denom = "numeric",
                                               w = "bigq", w_dyad = "bigq", approx_error = "numeric", tree = "Rdict",
                                               cone_lb = "numeric", cone_num = "numeric", cone_num_over = "numeric"),
                                prototype(p = NA_real_, max_denom = 1024), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @param max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @rdname GeoMean-class
GeoMean <- function(x, p = NA_real_, max_denom = 1024) { .GeoMean(x = x, p = p, max_denom  = max_denom) }

setMethod("initialize", "GeoMean", function(.Object, ..., x, p, max_denom) {
  .Object@x <- x
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x))

  x <- .Object@args[[1]]
  if(is_vector(x))
    n <- ifelse(1, ndim(x) == 0, max(dim(x)))
  else
    stop("x must be a row or column vector.")
  
  if(is.na(p))
    p <- rep(1, n)
  .Object@p <- p
  
  if(length(p) != n)
    stop("x and p must have the same number of elements.")

  if(any(p < 0) || sum(p) <= 0)
    stop("powers must be nonnegative and not all zero.")

  frac <- fracify(p, max_denom)
  .Object@w <- frac[[1]]
  .Object@w_dyad <- frac[[2]]
  .Object@approx_error <- approx_error(p, .Object@w)

  .Object@tree <- decompose(.Object@w_dyad)

  # known lower bound on number of cones needed to represent w_dyad
  .Object@cone_lb <- lower_bound(.Object@w_dyad)

  # number of cones used past known lower bound
  .Object@cone_num_over <- over_bound(.Object@w_dyad, .Object@tree)

  # number of cones used
  .Object@cone_num <- .Object@cone_lb + .Object@cone_num_over
  .Object
})

#' @param object A \linkS4class{GeoMean} object.
#' @param values A list of arguments to the atom.
#' @describeIn GeoMean The (weighted) geometric mean of the elements of \code{x}.
setMethod("to_numeric", "GeoMean", function(object, values) {
  values <- as.numeric(values[[1]])
  val <- 1.0
  for(idx in 1:length(values)) {
    x <- values[[idx]]
    p <- object@w[idx]
    val <- val * Rmpfr::mpfr(x, Rmpfr::getPrec(x))^p
  }
  return(gmp::asNumeric(val))   # TODO: Handle mpfr objects in the backend later
})

setMethod(".domain", "GeoMean", function(object) { list(object@args[[1]][object@w > 0] >= 0) })

setMethod(".grad", "GeoMean", function(object, values) {
  x <- as.matrix(values[[1]])
  # No special case when only one non-zero weight
  w_arr <- as.numeric(object@w)
  # Outside domain
  if(any(x[w_arr > 0] <= 0))
    return(list(NA_real_))
  else {
    D <- w_arr/as.numeric(x) * to_numeric(object, values)
    return(list(Matrix(D, sparse = TRUE)))
  }
})

#' @describeIn GeoMean The atom is a scalar.
setMethod("dim_from_args", "GeoMean", function(object) { c() })

#' @describeIn GeoMean The atom is non-negative.
setMethod("sign_from_args", "GeoMean", function(object) { c(TRUE, FALSE) })

#' @describeIn GeoMean The atom is not convex.
setMethod("is_atom_convex", "GeoMean", function(object) { FALSE })

#' @describeIn GeoMean The atom is concave.
setMethod("is_atom_concave", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "GeoMean", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn GeoMean The atom is weakly increasing in every argument.
setMethod("is_incr", "GeoMean", function(object, idx) { TRUE })

#' @describeIn GeoMean The atom is not weakly decreasing in any argument.
setMethod("is_decr", "GeoMean", function(object, idx) { FALSE })

#' @describeIn GeoMean Returns \code{list(w, dyadic completion, tree of dyads)}.
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

setMethod("copy", "GeoMean", function(object, args = NULL, id_objects = list()) {
  if(is.null(args))
    args <- object@args
  
  copy <- do.call(class(object), args)
  data <- get_data(object)
  copy@w <- data[[1]]
  copy@w_dyad <- data[[2]]
  copy@tree <- data[[3]]
  
  copy@approx_error <- object@approx_error
  copy@cone_lb <- object@cone_lb
  copy@cone_num_over <- object@cone_num_over
  copy@cone_num <- object@cone_num
  copy
})

HarmonicMean <- function(x) {
  x <- as.Constant(x)
  size(x) * Pnorm(x = x, p = -1)
}

#'
#' The LambdaMax class.
#'
#' The maximum eigenvalue of a matrix, \eqn{\lambda_{\max}(A)}.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name LambdaMax-class
#' @aliases LambdaMax
#' @rdname LambdaMax-class
.LambdaMax <- setClass("LambdaMax", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname LambdaMax-class
LambdaMax <- function(A) { .LambdaMax(A = A) }

setMethod("initialize", "LambdaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{LambdaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn LambdaMax The largest eigenvalue of \code{A}. Requires that \code{A} be symmetric.
setMethod("to_numeric", "LambdaMax", function(object, values) {
  # if(any(t(values[[1]]) != values[[1]]))
  #  stop("LambdaMax called on a non-symmetric matrix")
  max(eigen(values[[1]], only.values = TRUE)$values)
})

#' @rdname LambdaMax Returns the constraints describing the domain of the atom.
setMethod(".domain", "LambdaMax", function(object) { list(Conj(t(object@args[[1]])) == object@args[[1]]) })

#' @rdname LambdaMax Gives the (sub/super)gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix.
setMethod(".grad", "LambdaMax", function(object, values) {
  r <- eigen(values[[1]], only.values = FALSE)
  v <- r$vectors  # eigenvectors
  w <- r$values   # eigenvalues
  
  d <- array(0, dim = dim(w))
  d[length(d)] <- 1
  d <- diag(d)
  D <- v %*% d %*% t(v)
  list(Matrix(as.numeric(D), sparse = TRUE))
})

#' @describeIn LambdaMax Check that \code{A} is square.
setMethod("validate_args", "LambdaMax", function(object) {
  if(ndim(object@args[[1]]) != 2 || nrow(object@args[[1]]) != ncol(object@args[[1]]))
    stop("The argument to LambdaMax must resolve to a square matrix")
})

#' @describeIn LambdaMax The atom is a scalar.
setMethod("dim_from_args", "LambdaMax", function(object) { c() })

#' @describeIn LambdaMax The sign of the atom is unknown.
setMethod("sign_from_args", "LambdaMax", function(object) { c(FALSE, FALSE) })

#' @describeIn LambdaMax The atom is convex.
setMethod("is_atom_convex", "LambdaMax", function(object) { TRUE })

#' @describeIn LambdaMax The atom is not concave.
setMethod("is_atom_concave", "LambdaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LambdaMax The atom is not monotonic in any argument.
setMethod("is_incr", "LambdaMax", function(object, idx) { FALSE })

#' @describeIn LambdaMax The atom is not monotonic in any argument.
setMethod("is_decr", "LambdaMax", function(object, idx) { FALSE })

LambdaMin <- function(A) {
  A <- as.Constant(A)
  -LambdaMax(-A)
}

#'
#' The LambdaSumLargest class.
#'
#' This class represents the sum of the \code{k} largest eigenvalues of a matrix.
#'
#' @slot k A positive integer.
#' @name LambdaSumLargest-class
#' @aliases LambdaSumLargest
#' @rdname LambdaSumLargest-class
.LambdaSumLargest <- setClass("LambdaSumLargest", representation(k = "numeric"), contains = "LambdaMax")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param k A positive integer.
#' @rdname LambdaSumLargest-class
LambdaSumLargest <- function(A, k) { .LambdaSumLargest(A = A, k = k) }

setMethod("initialize", "LambdaSumLargest", function(.Object, ..., k) {
  .Object@k <- k
  callNextMethod(.Object, ...)
})

#' @describeIn LambdaSumLargest Does the atom handle complex numbers?
setMethod("allow_complex", "LambdaSumLargest", function(object) { TRUE })

#' @param object A \linkS4class{LambdaSumLargest} object.
#' @param values A list of arguments to the atom.
#' @rdname LambdaSumLargest Returns the largest eigenvalue of \code{A}, which must be symmetric.
setMethod("to_numeric", "LambdaSumLargest", function(object, values) {
  # if(any(t(values[[1]]) != values[[1]]))
  #  stop("LambdaSumLargest called on a non-symmetric matrix")
  eigs <- eigen(values[[1]], only.values = TRUE)$values
  value(SumLargest(eigs, object@k))
})

#' @rdname LambdaSumLargest Verify that the argument \code{A} is square.
setMethod("validate_args", "LambdaSumLargest", function(object) {
  A <- object@args[[1]]
  if(ndim(A) != 2 || nrow(A) != ncol(A))
    stop("First argument must be a square matrix.")
  else if(as.integer(object@k) != object@k || object@k <= 0)
    stop("Second argument must be a positive integer.")
})

#' @rdname LambdaSumLargest Returns the parameter \code{k}.
setMethod("get_data", "LambdaSumLargest", function(object) { list(object@k) })
setMethod(".grad", "LambdaSumLargest", function(object, values) { stop("Unimplemented") })

LambdaSumSmallest <- function(A, k) {
  A <- as.Constant(A)
  -LambdaSumLargest(-A, k)
}

#'
#' The LogDet class.
#'
#' The natural logarithm of the determinant of a matrix, \eqn{\log\det(A)}.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name LogDet-class
#' @aliases LogDet
#' @rdname LogDet-class
.LogDet <- setClass("LogDet", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname LogDet-class
LogDet <- function(A) { .LogDet(A = A) }

setMethod("initialize", "LogDet", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{LogDet} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogDet The log-determinant of SDP matrix \code{A}. This is the sum of logs of the eigenvalues and is equivalent to the nuclear norm of the matrix logarithm of \code{A}.
setMethod("to_numeric", "LogDet", function(object, values) {
  logdet <- determinant(values[[1]], logarithm = TRUE)
  if(logdet$sign == 1)
    return(as.numeric(logdet$modulus))
  else
    return(-Inf)
})

#' @describeIn LogDet Check that \code{A} is square.
setMethod("validate_args", "LogDet", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(length(arg_dim) == 1 || arg_dim[1] != arg_dim[2])
    stop("The argument to LogDet must be a square matrix")
})

#' @describeIn LogDet The atom is a scalar.
setMethod("dim_from_args", "LogDet", function(object) { c() })

#' @describeIn LogDet The atom is non-negative.
setMethod("sign_from_args",  "LogDet", function(object) { c(TRUE, FALSE) })

#' @describeIn LogDet The atom is not convex.
setMethod("is_atom_convex", "LogDet", function(object) { FALSE })

#' @describeIn LogDet The atom is concave.
setMethod("is_atom_concave", "LogDet", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_incr", "LogDet", function(object, idx) { FALSE })

#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_decr", "LogDet", function(object, idx) { FALSE })

setMethod(".grad", "LogDet", function(object, values) {
  X <- as.matrix(values[[1]])
  eigen_val <- eigen(X, only.values = TRUE)$values
  if(min(eigen_val) > 0) {
    # Grad: t(X^(-1))
    D <- t(base::solve(X))
    return(list(Matrix(as.numeric(D), sparse = TRUE)))
  } else   # Outside domain
    return(list(NA_real_))
})

setMethod(".domain", "LogDet", function(object) { list(object@args[[1]] %>>% 0) })

#'
#' The LogSumExp class.
#'
#' The natural logarithm of the sum of the elementwise exponential, \eqn{\log\sum_{i=1}^n e^{x_i}}.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name LogSumExp-class
#' @aliases LogSumExp
#' @rdname LogSumExp-class
.LogSumExp <- setClass("LogSumExp", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname LogSumExp-class
LogSumExp <- function(x, axis = NA_real_) { .LogSumExp(expr = x, axis = axis) }

#' @param object A \linkS4class{LogSumExp} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogSumExp Evaluates \eqn{e^x} elementwise, sums, and takes the natural log.
setMethod("to_numeric", "LogSumExp", function(object, values) {
  if(is.na(object@axis))
    log(sum(exp(values[[1]])))
  else
    log(apply(exp(values[[1]]), object@axis, sum))
})

setMethod(".grad", "LogSumExp", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "LogSumExp", function(object, value) {
  denom <- sum(exp(value))
  nom <- exp(value)
  D <- nom/denom
  D
})

#' @describeIn LogSumExp Returns sign (is positive, is negative) of the atom.
setMethod("sign_from_args",  "LogSumExp", function(object) { c(FALSE, FALSE) })

#' @describeIn LogSumExp The atom is convex.
setMethod("is_atom_convex", "LogSumExp", function(object) { TRUE })

#' @describeIn LogSumExp The atom is not concave.
setMethod("is_atom_concave", "LogSumExp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is weakly increasing in the index.
setMethod("is_incr", "LogSumExp", function(object, idx) { TRUE })

#' @describeIn LogSumExp The atom is not weakly decreasing in the index.
setMethod("is_decr", "LogSumExp", function(object, idx) { FALSE })

#'
#' The MatrixFrac class.
#'
#' The matrix fraction function \eqn{tr(X^T P^{-1} X)}.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @slot P An \linkS4class{Expression} or numeric matrix.
#' @name MatrixFrac-class
#' @aliases MatrixFrac
#' @rdname MatrixFrac-class
.MatrixFrac <- setClass("MatrixFrac", representation(X = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @param P An \linkS4class{Expression} or numeric matrix.
#' @rdname MatrixFrac-class
MatrixFrac <- function(X, P) {
  if(is.array(P) || is.matrix(P)) {
    invP <- as.matrix(base::solve(P))
    .QuadForm(x = X, P = (invP + t(Conj(invP))) / 2.0)
  } else
    .MatrixFrac(X = X, P = P)
}

matrix_frac <- function(X, P) {
  if(is.array(P)) {
    invP <- base::solve(P)
    return(QuadForm(X, (invP + t(Conj(invP)))/2.0))
  } else
    return(MatrixFrac(X, P))
}

setMethod("initialize", "MatrixFrac", function(.Object, ..., X, P) {
  .Object@X <- X
  .Object@P <- P
  callNextMethod(.Object, ..., atom_args = list(.Object@X, .Object@P))
})

#' @describeIn MatrixFrac Does the atom handle complex numbers?
setMethod("allow_complex", "MatrixFrac", function(object) { TRUE })

#' @param object A \linkS4class{MatrixFrac} object.
#' @param values A list of arguments to the atom.
#' @describeIn MatrixFrac The trace of \eqn{X^TP^{-1}X}.
setMethod("to_numeric", "MatrixFrac", function(object, values) {
  # TODO: Raise error if not invertible?
  X <- values[[1]]
  P <- values[[2]]
  if(is_complex(object@args[[1]]))
    product <- t(Conj(X)) %*% base::solve(P) %*% X
  else
    product <- t(X) %*% base::solve(P) %*% X
  
  if(length(dim(product)) == 2)
    return(sum(diag(product)))
  else
    return(product)
})

#' @describeIn MatrixFrac Check that the dimensions of \code{x} and \code{P} match.
setMethod("validate_args", "MatrixFrac", function(object) {
  X <- object@args[[1]]
  P <- object@args[[2]]
  if(ndim(P) != 2 || nrow(P) != ncol(P))
    stop("The second argument to MatrixFrac must be a square matrix.")
  else if(nrow(X) != nrow(P))
    stop("The arguments to MatrixFrac have incompatible dimensions.")
})

#' @describeIn MatrixFrac The atom is a scalar.
setMethod("dim_from_args", "MatrixFrac", function(object) { c() })

#' @describeIn MatrixFrac The atom is positive.
setMethod("sign_from_args", "MatrixFrac", function(object) { c(TRUE, FALSE) })

#' @describeIn MatrixFrac The atom is convex.
setMethod("is_atom_convex", "MatrixFrac", function(object) { TRUE })

#' @describeIn MatrixFrac The atom is not concave.
setMethod("is_atom_concave", "MatrixFrac", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MatrixFrac The atom is not monotonic in any argument.
setMethod("is_incr", "MatrixFrac", function(object, idx) { FALSE })

#' @describeIn MatrixFrac The atom is not monotonic in any argument.
setMethod("is_decr", "MatrixFrac", function(object, idx) { FALSE })

#' @describeIn MatrixFrac True if x is affine and P is constant.
setMethod("is_quadratic", "MatrixFrac", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn MatrixFrac True if x is piecewise linear and P is constant.
setMethod("is_qpwa", "MatrixFrac", function(object) { is_pwl(object@args[[1]]) && is_constant(object@args[[2]]) })

setMethod(".domain", "MatrixFrac", function(object) { list(object@args[[2]] %>>% 0) })

setMethod(".grad", "MatrixFrac", function(object, values) {
  X <- as.matrix(values[[1]])
  P <- as.matrix(values[[2]])
  P_inv <- tryCatch({
    base::solve(P)
  }, error = function(e) {
    return(list(NA_real_, NA_real_))
  })

  if(is.null(dim(P_inv)) && is.na(P_inv))
    return(list(NA_real_, NA_real_))

  # partial_X = (P^-1+P^-T)X
  # partial_P = (P^-1 * X * X^T * P^-1)^T
  DX <- (P_inv + t(P_inv)) %*% X
  DX <- as.numeric(t(DX))
  DX <- Matrix(DX, sparse = TRUE)

  DP <- P_inv %*% X
  DP <- DP %*% t(X)
  DP <- DP %*% P_inv
  DP <- -t(DP)
  DP <- Matrix(as.numeric(t(DP)), sparse = TRUE)
  list(DX, DP)
})

#'
#' The MaxEntries class.
#'
#' The maximum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name MaxEntries-class
#' @aliases MaxEntries
#' @rdname MaxEntries-class
.MaxEntries <- setClass("MaxEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname MaxEntries-class
MaxEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MaxEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MaxEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MaxEntries The largest entry in \code{x}.
setMethod("to_numeric", "MaxEntries", function(object, values) {
  apply_with_keepdims(values[[1]], max, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MaxEntries The sign of the atom.
setMethod("sign_from_args",  "MaxEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MaxEntries The atom is convex.
setMethod("is_atom_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries The atom is not concave.
setMethod("is_atom_concave", "MaxEntries", function(object) { FALSE })

#' @describeIn MaxEntries Is the atom log-log convex.
setMethod("is_atom_log_log_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries Is the atom log-log concave.
setMethod("is_atom_log_log_concave", "MaxEntries", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MaxEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MaxEntries", function(object, idx) { TRUE })

#' @describeIn MaxEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MaxEntries", function(object, idx) { FALSE })

#' @describeIn MaxEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MaxEntries", function(object) { is_pwl(object@args[[1]]) })

setMethod(".grad", "MaxEntries", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "MaxEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.numeric(value)
  idx <- (value == max(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

#'
#' The MinEntries class.
#'
#' The minimum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name MinEntries-class
#' @aliases MinEntries
#' @rdname MinEntries-class
.MinEntries <- setClass("MinEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname MinEntries-class
MinEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MinEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MinEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MinEntries The largest entry in \code{x}.
setMethod("to_numeric", "MinEntries", function(object, values) {
  apply_with_keepdims(values[[1]], min, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MinEntries The sign of the atom.
setMethod("sign_from_args",  "MinEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MinEntries The atom is not convex.
setMethod("is_atom_convex", "MaxEntries", function(object) { FALSE })

#' @describeIn MinEntries The atom is concave.
setMethod("is_atom_concave", "MinEntries", function(object) { TRUE })

#' @describeIn MinEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MinEntries", function(object) { FALSE })

#' @describeIn MinEntries Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MinEntries", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MinEntries", function(object, idx) { TRUE })

#' @describeIn MinEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MinEntries", function(object, idx) { FALSE })

#' @describeIn MinEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MinEntries", function(object) { is_pwl(object@args[[1]]) })

setMethod(".grad", "MinEntries", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "MinEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.numeric(value)
  idx <- (value == min(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

#'
#' The Pnorm class.
#'
#' This class represents the vector p-norm.
#'
#' If given a matrix variable, \code{Pnorm} will treat it as a vector and compute the p-norm of the concatenated columns.
#'
#' For \eqn{p \geq 1}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n}.
#' For \eqn{p < 1, p\neq 0}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n x_i^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n_+}.
#'
#' \itemize{
#'    \item Note that the "p-norm" is actually a \strong{norm} only when \eqn{p \geq 1} or \eqn{p = +\infty}. For these cases, it is convex.
#'    \item The expression is undefined when \eqn{p = 0}.
#'    \item Otherwise, when \eqn{p < 1}, the expression is concave, but not a true norm.
#' }
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot p A number greater than or equal to 1, or equal to positive infinity.
#' @slot max_denom The maximum denominator considered in forming a rational approximation for \eqn{p}.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot .approx_error (Internal) The absolute difference between \eqn{p} and its rational approximation.
#' @name Pnorm-class
#' @aliases Pnorm
#' @rdname Pnorm-class
.Pnorm <- setClass("Pnorm", representation(p = "numeric", max_denom = "numeric", .approx_error = "numeric", .original_p = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_, .original_p = NA_real_), contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param p A number greater than or equal to 1, or equal to positive infinity.
#' @param max_denom The maximum denominator considered in forming a rational approximation for \eqn{p}.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname Pnorm-class
Pnorm <- function(x, p = 2, axis = NA_real_, keepdims = FALSE, max_denom = 1024) {
  if(p == 1)
    Norm1(x, axis = axis, keepdims = keepdims)
  else if(p %in% c(Inf, "inf", "Inf"))
    NormInf(x, axis = axis, keepdims = keepdims)
  else
    .Pnorm(expr = x, axis = axis, keepdims = keepdims, p = p, max_denom = max_denom)
}

setMethod("initialize", "Pnorm", function(.Object, ..., p = 2, max_denom = 1024, .approx_error = NA_real_, .original_p = NA_real_) {
  if(p == 1)
    stop("Use the Norm1 class to instantiate a 1-norm.")
  else if(p %in% c(Inf, "inf", "Inf"))
    stop("Use the NormInf class to instantiate an infinity-norm.")
  # else if(p < 0)
  #  .Object@p <- pow_neg(p, max_denom)
  # else if(p > 0 && p < 1)
  #  .Object@p <- pow_mid(p, max_denom)
  # else if(p > 1)
  #  .Object@p <- pow_high(p, max_denom)
  # else
  #  stop("Invalid value of p.")
  .Object@p <- p
  
  .Object@max_denom <- max_denom
  .Object@.approx_error <- abs(.Object@p - p)
  .Object@.original_p <- p
  callNextMethod(.Object, ...)
})

# Internal method for calculating the p-norm
.p_norm <- function(x, p) {
  if(p == Inf)
    max(abs(x))
  else if(p == 0)
    sum(x != 0)
  else if(p %% 2 == 0 || p < 1)
    sum(x^p)^(1/p)
  else if(p >= 1)
    sum(abs(x)^p)^(1/p)
  else
    stop("Invalid p = ", p)
}

#' @describeIn Pnorm Does the atom handle complex numbers?
setMethod("allow_complex", "Pnorm", function(object) { TRUE })

#' @param object A \linkS4class{Pnorm} object.
#' @param values A list of arguments to the atom.
#' @describeIn Pnorm The p-norm of \code{x}.
setMethod("to_numeric", "Pnorm", function(object, values) {
  if(is.na(object@axis))
    values <- as.numeric(values[[1]])
  else
    values <- as.matrix(values[[1]])
  
  if(object@p < 1 && any(values < 0))
    return(-Inf)
  if(object@p < 0 && any(values == 0))
    return(0)
  
  apply_with_keepdims(values, function(x) { .p_norm(x, object@p) }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn Pnorm Check that the arguments are valid.
setMethod("validate_args", "Pnorm", function(object) {
  callNextMethod()
  if(!is.na(object@axis) && object@p != 2)
    stop("The axis parameter is only supported for p = 2.")
  if(object@p < 1 && is_complex(object@args[[1]]))
    stop("Pnorm(x, p) cannot have x complex for p < 1.")
})

#' @describeIn Pnorm The atom is positive.
setMethod("sign_from_args",  "Pnorm", function(object) { c(TRUE, FALSE) })

#' @describeIn Pnorm The atom is convex if \eqn{p \geq 1}.
setMethod("is_atom_convex", "Pnorm", function(object) { object@p > 1 })

#' @describeIn Pnorm The atom is concave if \eqn{p < 1}.
setMethod("is_atom_concave", "Pnorm", function(object) { object@p < 1 })

#' @describeIn Pnorm Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Pnorm", function(object) { TRUE })

#' @describeIn Pnorm Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Pnorm", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Pnorm The atom is weakly increasing if \eqn{p < 1} or \eqn{p > 1} and \code{x} is positive.
setMethod("is_incr", "Pnorm", function(object, idx) { object@p < 1 || (object@p > 1 && is_nonneg(object@args[[1]])) })

#' @describeIn Pnorm The atom is weakly decreasing if \eqn{p > 1} and \code{x} is negative.
setMethod("is_decr", "Pnorm", function(object, idx) { object@p > 1 && is_nonpos(object@args[[1]]) })

#' @describeIn Pnorm The atom is not piecewise linear unless \eqn{p = 1} or \eqn{p = \infty}.
setMethod("is_pwl", "Pnorm", function(object) { FALSE })

#' @describeIn Pnorm Returns \code{list(p, axis)}.
setMethod("get_data", "Pnorm", function(object) { list(object@p, object@axis) })

#' @describeIn Pnorm The name and arguments of the atom.
setMethod("name", "Pnorm", function(x) {
  sprintf("%s(%s, %s)", class(x), name(x@args[[1]]), x@p)
})

setMethod(".domain", "Pnorm", function(object) {
  if(object@p < 1 && object@p != 0)
    list(object@args[[1]] >= 0)
  else
    list()
})

setMethod(".grad", "Pnorm", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "Pnorm", function(object, value) {
  rows <- size(object@args[[1]])
  value <- as.matrix(value)

  # Outside domain
  if(object@p < 1 && any(value <= 0))
    return(NA_real_)
  
  D_null <- sparseMatrix(i = c(), j = c(), dims = c(rows, 1))
  denominator <- .p_norm(value, object@p)
  denominator <- denominator^(object@p - 1)

  # Subgrad is 0 when denom is 0 (or undefined)
  if(denominator == 0) {
    if(object@p > 1)
      return(D_null)
    else
      return(NA_real_)
  } else {
    numerator <- value^(object@p - 1)
    frac <- numerator / denominator
    return(matrix(as.numeric(frac)))
  }
})

MixedNorm <- function(X, p = 2, q = 1) {
  X <- as.Constant(X)
  
  # Inner norms
  vecnorms <- Norm(X, p, axis = 1)
  
  # Outer norms
  Norm(vecnorms, q)
}

Norm <- function(x, p = 2, axis = NA_real_) {
  x <- as.Constant(x)
  
  # Matrix norms take precedence.
  num_nontrivial_idxs <- sum(dim(x) > 1)
  if(is.na(axis) && ndim(x) == 2) {
    if(p == 1)   # Matrix 1-norm.
      MaxEntries(Norm1(x, axis = 2))
    else if(p == "fro" || (p == 2 && num_nontrivial_idxs == 1))   # Frobenius norm.
      Pnorm(Vec(x), 2)
    else if(p == 2)   # Matrix 2-norm is largest singular value.
      SigmaMax(x)
    else if(p == "nuc")   # The nuclear norm (sum of singular values)
      NormNuc(x)
    else if(p %in% c(Inf, "inf", "Inf"))   # The matrix infinity-norm.
      MaxEntries(Norm1(x, axis = 1))
    else
      stop("Unsupported matrix norm.")
  } else {
    if(p == 1 || is_scalar(x))
      Norm1(x, axis = axis)
    else if(p %in% c(Inf, "inf", "Inf"))
      NormInf(x, axis)
    else
      Pnorm(x, p, axis)
  }
}

Norm2 <- function(x, axis = NA_real_, keepdims = FALSE, max_denom = 1024) {
  Pnorm(x, p = 2, axis = axis, keepdims = keepdims, max_denom = max_denom)
}

#'
#' The Norm1 class.
#' 
#' This class represents the 1-norm of an expression.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Norm1-class
#' @aliases Norm1
#' @rdname Norm1-class
.Norm1 <- setClass("Norm1", contains = "AxisAtom")

Norm1 <- function(x, axis = NA_real_, keepdims = FALSE) { .Norm1(x = x, axis = axis, keepdims = keepdims) }

#' @param x,object A \linkS4class{Norm1} object.
setMethod("name", "Norm1", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ")", sep = "")
})

#' @param values A list of arguments to the atom.
#' @rdname Norm1 Returns the 1-norm of x along the given axis.
setMethod("to_numeric", "Norm1", function(object, values) {
  if(is.na(object@axis))
    base::norm(values[[1]], type = "O")
  else
    apply_with_keepdims(values[[1]], function(x) { norm(x, type = "O") }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn Norm1 Does the atom handle complex numbers?
setMethod("allow_complex", "Norm1", function(object) { TRUE })

#' @rdname Norm1 The atom is always positive.
setMethod("sign_from_args", "Norm1", function(object) { c(TRUE, FALSE) })

#' @rdname Norm1 The atom is convex.
setMethod("is_atom_convex", "Norm1", function(object) { TRUE })

#' @rdname Norm1 The atom is not concave.
setMethod("is_atom_concave", "Norm1", function(object) { FALSE })

#' @rdname Norm1 Is the composition weakly increasing in argument \code{idx}?
setMethod("is_incr", "Norm1", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @rdname Norm1 Is the composition weakly decreasing in argument \code{idx}?
setMethod("is_decr", "Norm1", function(object, idx) { is_nonpos(object@args[[1]]) })

#' @rdname Norm1 Is the atom piecewise linear?
setMethod("is_pwl", "Norm1", function(object) {
  is_pwl(object@args[[1]]) && (is_real(object@args[[1]]) || is_imag(object@args[[1]]))
})

#' @rdname Norm1 Returns the axis.
setMethod("get_data", "Norm1", function(object) { list(object@axis) })

setMethod(".domain", "Norm1", function(object) { list() })
setMethod(".grad", "Norm1", function(object, values) { .axis_grad(object, values) })
setMethod(".column_grad", "Norm1", function(object, value) {
  rows <- size(object@args[[1]])
  D_null <- Matrix(0, nrow = rows, ncol = 1, sparse = TRUE)
  D_null <- D_null + (value > 0)
  D_null <- D_null - (value < 0)
  D_null   # TODO: Check this is same as ravel and transpose command in CVXPY.
})

#'
#' The NormInf class.
#'
#' This class represents the infinity-norm.
#'
#' @name NormInf-class
#' @aliases NormInf
#' @rdname NormInf-class
NormInf <- setClass("NormInf", contains = "AxisAtom")

setMethod("name", "NormInf", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ")", sep = "")
})

#' @describeIn NormInf Returns the infinity norm of \code{x}.
setMethod("to_numeric", "NormInf", function(object, values) {
  if(is.na(object@axis))
    base::norm(values[[1]], type = "I")
  else
    apply_with_keepdims(values[[1]], function(x) { norm(x, type = "I") }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn NormInf Does the atom handle complex numbers?
setMethod("allow_complex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf The atom is always positive.
setMethod("sign_from_args", "NormInf", function(object) { c(TRUE, FALSE) })

#' @describeIn NormInf The atom is convex.
setMethod("is_atom_convex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf The atom is not concave.
setMethod("is_atom_concave", "NormInf", function(object) { FALSE })

#' @describeIn NormInf Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "NormInf", function(object) { FALSE })

#' @describeIn NormInf Is the composition weakly increasing in argument \code{idx}?
setMethod("is_incr", "NormInf", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @describeIn NormInf Is the composition weakly decreasing in argument \code{idx}?
setMethod("is_decr", "NormInf", function(object, idx) { is_nonpos(object@args[[1]]) })

#' @describeIn NormInf Is the atom piecewise linear?
setMethod("is_pwl", "NormInf", function(object) { is_pwl(object@args[[1]]) })

#' @describeIn NormInf Returns the axis.
setMethod("get_data", "NormInf", function(object) { list(object@axis) })

setMethod(".domain", "NormInf", function(object) { list() })
setMethod(".grad", "NormInf", function(object, values) { .axis_grad(object, values) })
setMethod(".column_grad", "NormInf", function(object, value) { stop("Unimplemented") })   # TODO: Implement this! })

#'
#' The NormNuc class.
#'
#' The nuclear norm, i.e. sum of the singular values of a matrix.
#'
#' @slot A An \linkS4class{Expression} representing a matrix.
#' @name NormNuc-class
#' @aliases NormNuc
#' @rdname NormNuc-class
.NormNuc <- setClass("NormNuc", representation(A = "Expression"), contains = "Atom")

#' @param A An \linkS4class{Expression} representing a matrix.
#' @rdname NormNuc-class
NormNuc <- function(A) { .NormNuc(A = A) }

setMethod("initialize", "NormNuc", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{NormNuc} object.
#' @param values A list of arguments to the atom.
#' @describeIn NormNuc The nuclear norm (i.e., the sum of the singular values) of \code{A}.
setMethod("to_numeric", "NormNuc", function(object, values) {
  # Returns the nuclear norm (i.e. the sum of the singular values) of A
  sum(svd(values[[1]])$d)
})

#' @describeIn NormNuc Does the atom handle complex numbers?
setMethod("allow_complex", "NormNuc", function(object) { TRUE })

#' @describeIn NormNuc The atom is a scalar.
setMethod("dim_from_args", "NormNuc", function(object) { c() })

#' @describeIn NormNuc The atom is positive.
setMethod("sign_from_args",  "NormNuc", function(object) { c(TRUE, FALSE) })

#' @describeIn NormNuc The atom is convex.
setMethod("is_atom_convex", "NormNuc", function(object) { TRUE })

#' @describeIn NormNuc The atom is not concave.
setMethod("is_atom_concave", "NormNuc", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn NormNuc The atom is not monotonic in any argument.
setMethod("is_incr", "NormNuc", function(object, idx) { FALSE })

#' @describeIn NormNuc The atom is not monotonic in any argument.
setMethod("is_decr", "NormNuc", function(object, idx) { FALSE })

setMethod(".grad", "NormNuc", function(object, values) {
  # Grad: UV^T
  s <- svd(values[[1]])
  D <- s$u %*% t(s$v)
  list(Matrix(as.numeric(D), sparse = TRUE))
})

#'
#' The OneMinusPos class.
#'
#' This class represents the difference \eqn{1 - x} with domain \eqn{\{x : 0 < x < 1}\}
#' 
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @name OneMinusPos-class
#' @aliases OneMinusPos
#' @rdname OneMinusPos-class
.OneMinusPos <- setClass("OneMinusPos", representation(x = "ConstValORExpr", .ones = "numeric"), prototype(.ones = NA_real_), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @rdname OneMinusPos-class
OneMinusPos <- function(x) { .OneMinusPos(x = x) }

setMethod("initialize", "OneMinusPos", function(.Object, ..., x) {
  .Object@x <- x
  .Object@.ones <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x))
  .Object@args[[1]] <- x
  .Object
})

#' @param x,object A \linkS4class{OneMinusPos} object.
setMethod("name", "OneMinusPos", function(x) { paste(class(x), x@args[[1]]) })

#' @param values A list of arguments to the atom.
#' @describeIn OneMinusPos Returns one minus the value.
setMethod("to_numeric", "OneMinusPos", function(object, values) { object@.ones - values[[1]] })

#' @describeIn OneMinusPos The dimensions of the atom.
setMethod("dim_from_args", "OneMinusPos", function(object) { dim(object@args[[1]]) })

#' @describeIn OneMinusPos Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "OneMinusPos", function(object) { c(TRUE, FALSE) })

#' @describeIn OneMinusPos Is the atom convex?
setMethod("is_atom_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom concave?
setMethod("is_atom_concave", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "OneMinusPos", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn OneMinusPos Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "OneMinusPos", function(object, idx) { FALSE })

#' @describeIn OneMinusPos Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "OneMinusPos", function(object, idx) { TRUE })

setMethod(".grad", "OneMinusPos", function(object, values) { Matrix(-object@.ones, sparse = TRUE) })

# The difference x - y with domain {x,y: x > y > 0}.
DiffPos <- function(x, y) {
  x * OneMinusPos(y/x)
}

#'
#' The PfEigenvalue class.
#'
#' This class represents the Perron-Frobenius eigenvalue of a positive matrix.
#' 
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @name PfEigenvalue-class
#' @aliases PfEigenvalue
#' @rdname PfEigenvalue-class
.PfEigenvalue <- setClass("PfEigenvalue", representation(X = "ConstValORExpr"), 
                          validity = function(object) {
                            if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                              stop("[PfEigenvalue: X] The argument X must be a square matrix")
                            return(TRUE)
                          }, contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @rdname PfEigenvalue-class
PfEigenvalue <- function(X) { .PfEigenvalue(X = X) }

setMethod("initialize", "PfEigenvalue", function(.Object, ..., X = X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param x,object A \linkS4class{PfEigenvalue} object.
setMethod("name", "PfEigenvalue", function(x) { paste(class(x), x@args[[1]]) })

#' @param values A list of arguments to the atom.
#' @describeIn PfEigenvalue Returns the Perron-Frobenius eigenvalue of \code{X}.
setMethod("to_numeric", "PfEigenvalue", function(object, values) {
  eig <- eigen(values[[1]], only.values = TRUE)
  max(abs(eig$values))
})

#' @describeIn PfEigenvalue The dimensions of the atom.
setMethod("dim_from_args", "PfEigenvalue", function(object) { NULL })

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

#' @describeIn PfEigenvalue Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "PfEigenvalue", function(object, idx) { FALSE })

setMethod(".grad", "PfEigenvalue", function(object, values) { NA_real_ })

#'
#' The ProdEntries class.
#'
#' The product of the entries in an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name ProdEntries-class
#' @aliases ProdEntries
#' @rdname ProdEntries-class
.ProdEntries <- setClass("ProdEntries", contains = "AxisAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname ProdEntries-class
ProdEntries <- function(expr, axis = NA_real_, keepdims = FALSE) { .ProdEntries(expr = expr, axis = axis, keepdims = keepdims) }

#' @describeIn ProdEntries The product of all the entries.
setMethod("to_numeric", "ProdEntries", function(object, values) {
  apply_with_keepdims(values[[1]], prod, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn ProdEntries Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "ProdEntries", function(object) { 
  if(is_nonneg(object@args[[1]]))
    c(TRUE, FALSE)
  else
    c(FALSE, FALSE)
})

#' @describeIn ProdEntries Is the atom convex?
setMethod("is_atom_convex", "ProdEntries", function(object) { FALSE })

#' @describeIn ProdEntries Is the atom concave?
setMethod("is_atom_concave", "ProdEntries", function(object) { FALSE })

#' @describeIn ProdEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "ProdEntries", function(object) { TRUE })

#' @describeIn ProdEntries is the atom log-log concave?
setMethod("is_atom_log_log_concave", "ProdEntries", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn ProdEntries Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "ProdEntries", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @describeIn ProdEntries Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "ProdEntries", function(object, idx) { FALSE })

setMethod(".column_grad", "ProdEntries", function(object, value) { prod(value)/value })
setMethod(".grad", "ProdEntries", function(object, values) { .axis_grad(object, values) })

#'
#' The QuadForm class.
#'
#' This class represents the quadratic form \eqn{x^T P x}
#' 
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @name QuadForm-class
#' @aliases QuadForm
#' @rdname QuadForm-class
.QuadForm <- setClass("QuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @rdname QuadForm-class
QuadForm <- function(x, P) { .QuadForm(x = x, P = P) }

setMethod("initialize", "QuadForm", function(.Object, ..., x, P) {
  .Object@x <- x
  .Object@P <- P
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@P))
})

setMethod("name", "QuadForm", function(x) {
  paste(class(x), "(", x@args[[1]], ", ", x@args[[2]], ")", sep = "")
})

#' @describeIn QuadForm Does the atom handle complex numbers?
setMethod("allow_complex", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Returns the quadratic form.
setMethod("to_numeric", "QuadForm", function(object, values) {
  prod <- values[[2]] %*% values[[1]]
  if(is_complex(object@args[[1]]))
    return(t(Conj(values[[1]])) %*% prod)
  else
    return(t(values[[1]]) %*% prod)
})

#' @describeIn QuadForm Checks the dimensions of the arguments.
setMethod("validate_args", "QuadForm", function(object) {
  callNextMethod()
  n <- nrow(object@args[[2]])
  if(ncol(object@args[[2]]) != n || !(dim(object@args[[1]]) %in% list(c(n, 1), c(n, NA_real_))))
    stop("Invalid dimensions for arguments.")
})

#' @describeIn QuadForm Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "QuadForm", function(object) { c(is_atom_convex(object), is_atom_concave(object)) })

#' @describeIn QuadForm The dimensions of the atom.
setMethod("dim_from_args", "QuadForm", function(object) { 
  if(ndim(object@args[[1]]) == 0)
    c()
  else
    c(1,1)
})

#' @describeIn QuadForm Is the atom convex?
setMethod("is_atom_convex", "QuadForm", function(object) { is_psd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom concave?
setMethod("is_atom_concave", "QuadForm", function(object) { is_nsd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadForm", function(object) { FALSE })

#' @describeIn QuadForm Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonneg(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonneg(object@args[[2]]))
})

#' @describeIn QuadForm Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonpos(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonpos(object@args[[2]]))
})

#' @describeIn QuadForm Is the atom quadratic?
setMethod("is_quadratic", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom piecewise linear?
setMethod("is_pwl", "QuadForm", function(object) { FALSE })

setMethod(".grad", "QuadForm", function(object, values) {
  x <- values[[1]]
  P <- values[[2]]
  D <- 2*P %*% t(x)
  Matrix(as.vector(t(D)), sparse = TRUE)
})

#'
#' The SymbolicQuadForm class.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @slot original_expression The original \linkS4class{Expression}.
#' @name SymbolicQuadForm-class
#' @aliases SymbolicQuadForm
#' @rdname SymbolicQuadForm-class
.SymbolicQuadForm <- setClass("SymbolicQuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr", original_expression = "Expression"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @param original_expression The original \linkS4class{Expression}.
#' @rdname SymbolicQuadForm-class
SymbolicQuadForm <- function(x, P, expr) { .SymbolicQuadForm(x = x, P = P, original_expression = expr) }

setMethod("initialize", "SymbolicQuadForm", function(.Object, ..., x, P, original_expression) {
  .Object@original_expression <- original_expression
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@P))
  .Object@P <- .Object@args[[2]]
  .Object
})

#' @rdname SymbolicQuadForm The dimensions of the atom.
setMethod("dim_from_args", "SymbolicQuadForm", function(object) { dim_from_args(object@original_expression) })

#' @rdname SymbolicQuadForm The sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "SymbolicQuadForm", function(object) { sign_from_args(object@original_expression) })

#' @rdname SymbolicQuadForm The original expression.
setMethod("get_data", "SymbolicQuadForm", function(object) { list(object@original_expression) })

#' @rdname SymbolicQuadForm Is the original expression convex?
setMethod("is_atom_convex", "SymbolicQuadForm", function(object) { is_atom_convex(object@original_expression) })

#' @rdname SymbolicQuadForm Is the original expression concave?
setMethod("is_atom_concave", "SymbolicQuadForm", function(object) { is_atom_concave(object@original_expression) })

#' @rdname SymbolicQuadForm Is the original expression weakly increasing in argument \code{idx}?
setMethod("is_incr", "SymbolicQuadForm", function(object, idx) { is_incr(object@original_expression, idx) })

#' @rdname SymbolicQuadForm Is the original expression weakly decreasing in argument \code{idx}?
setMethod("is_decr", "SymbolicQuadForm", function(object, idx) { is_decr(object@original_expression, idx) })

#' @rdname SymbolicQuadForm The atom is quadratic.
setMethod("is_quadratic", "SymbolicQuadForm", function(object) { TRUE })

setMethod(".grad", "SymbolicQuadForm", function(object, values) { stop("Unimplemented") })

.decomp_quad <- function(P, cond = NA, rcond = NA) {
  eig <- eigen(P, only.values = FALSE)
  w <- eig$values
  V <- eig$vectors
  
  if(!is.na(rcond))
    cond <- rcond
  if(cond == -1 || is.na(cond))
    cond <- 1e6 * .Machine$double.eps   # All real numbers are stored as double precision in R
  
  scale <- max(abs(w))
  if(scale < cond)
    return(list(scale = 0, M1 = V[,FALSE], M2 = V[,FALSE]))
  w_scaled <- w / scale
  maskp <- w_scaled > cond
  maskn <- w_scaled < -cond
  
  # TODO: Allow indefinite QuadForm
  if(any(maskp) && any(maskn))
    warning("Forming a non-convex expression QuadForm(x, indefinite)")
  
  if(sum(maskp) <= 1)
    M1 <- as.matrix(V[,maskp] * sqrt(w_scaled[maskp]))
  else
    M1 <- V[,maskp] %*% diag(sqrt(w_scaled[maskp]))
  
  if(sum(maskn) <= 1)
    M2 <- as.matrix(V[,maskn]) * sqrt(-w_scaled[maskn])
  else
    M2 <- V[,maskn] %*% diag(sqrt(-w_scaled[maskn]))
  list(scale = scale, M1 = M1, M2 = M2)
}

QuadForm <- function(x, P) {
  # x^T P x
  x <- as.Constant(x)
  P <- as.Constant(P)
  
  # Check dimensions.
  P_dim <- dim(P)
  if(ndim(P) != 2 || P_dim[1] != P_dim[2] || max(nrow(x), 1) != P_dim[1])
    stop("Invalid dimensions for arguments.")
  
  # P cannot be a parameter.
  if(is_constant(x))
    Conj(t(x)) %*% P %*% x
  else if(is_constant(P))
    QuadForm(x, P)
  else
    stop("At least one argument to QuadForm must be constant.")
}

#'
#' The QuadOverLin class.
#'
#' This class represents the sum of squared entries in X divided by a scalar y, \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @slot y A scalar \linkS4class{Expression} or numeric constant.
#' @name QuadOverLin-class
#' @aliases QuadOverLin
#' @rdname QuadOverLin-class
.QuadOverLin <- setClass("QuadOverLin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param y A scalar \linkS4class{Expression} or numeric constant.
#' @rdname QuadOverLin-class
QuadOverLin <- function(x, y) { .QuadOverLin(x = x, y = y) }

setMethod("initialize", "QuadOverLin", function(.Object, ..., x = .Object@x, y = .Object@y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@y))
})

#' @param object A \linkS4class{QuadOverLin} object.
#' @param values A list of arguments to the atom.
#' @describeIn QuadOverLin The sum of the entries of \code{x} squared over \code{y}.
setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(values[[1]]^2) / values[[2]] })

#' @describeIn QuadOverLin Check the dimensions of the arguments.
setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@args[[2]]))
    stop("The second argument to QuadOverLin must be a scalar.")
  callNextMethod()
})

#' @describeIn QuadOverLin The atom is a scalar.
setMethod("dim_from_args", "QuadOverLin", function(object) { c() })

#' @describeIn QuadOverLin The atom is positive.
setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })

#' @describeIn QuadOverLin The atom is convex.
setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin The atom is not concave.
setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })

#' @describeIn QuadOverLin Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadOverLin", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_nonneg(object@args[[idx]]) })

#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_nonpos(object@args[[idx]])) || (idx == 2) })

#' @describeIn QuadOverLin Quadratic if \code{x} is affine and \code{y} is constant.
setMethod("is_quadratic", "QuadOverLin", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn QuadOverLin Quadratic of piecewise affine if \code{x} is piecewise linear and \code{y} is constant.
setMethod("is_qpwa", "QuadOverLin", function(object) { is_pwl(object@args[[1]]) && is_constant(object@args[[2]]) })

setMethod(".domain", "QuadOverLin", function(object) { list(object@args[[2]] >= 0) })

setMethod(".grad", "QuadOverLin", function(object, values) {
  X <- values[[1]]
  y <- values[[2]]
  if(y <= 0)
    return(list(NA_real_, NA_real_))
  else {
    # DX = 2X/y, Dy = -||X||^2_2/y^2
    Dy <- -sum(X^2)/y^2
    Dy <- Matrix(Dy, sparse = TRUE)
    DX <- 2.0*X/y
    DX <- Matrix(as.numeric(t(DX)), sparse = TRUE)
    return(list(DX, Dy))
  }
})

#'
#' The SigmaMax class.
#'
#' The maximum singular value of a matrix.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name SigmaMax-class
#' @aliases SigmaMax
#' @rdname SigmaMax-class
.SigmaMax <- setClass("SigmaMax", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or matrix.
#' @rdname SigmaMax-class
SigmaMax <- function(A = A) { .SigmaMax(A = A) }

setMethod("initialize", "SigmaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{SigmaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn SigmaMax The largest singular value of \code{A}.
setMethod("to_numeric", "SigmaMax", function(object, values) { base::norm(values[[1]], type = "2") })

#' @describeIn SigmaMax Does the atom handle complex numbers?
setMethod("allow_complex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is a scalar.
setMethod("dim_from_args", "SigmaMax", function(object) { c() })

#' @describeIn SigmaMax The atom is positive.
setMethod("sign_from_args",  "SigmaMax", function(object) { c(TRUE, FALSE) })

#' @describeIn SigmaMax The atom is convex.
setMethod("is_atom_convex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is concave.
setMethod("is_atom_concave", "SigmaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_incr", "SigmaMax", function(object, idx) { FALSE })

#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_decr", "SigmaMax", function(object, idx) { FALSE })

setMethod(".grad", "SigmaMax", function(object, values) {
  # Grad: U diag(e_1) t(V)
  s <- svd(values[[1]])
  ds <- rep(0, length(s$d))
  ds[1] <- 1
  D <- s$u %*% diag(ds) %*% t(s$v)
  list(Matrix(as.numeric(D), sparse = TRUE))
})

#'
#' The SumLargest class.
#'
#' The sum of the largest k values of a matrix.
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @slot k The number of largest values to sum over.
#' @name SumLargest-class
#' @aliases SumLargest
#' @rdname SumLargest-class
.SumLargest <- setClass("SumLargest", representation(x = "ConstValORExpr", k = "numeric"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param k The number of largest values to sum over.
#' @rdname SumLargest-class
SumLargest <- function(x, k) { .SumLargest(x = x, k = k) }

setMethod("initialize", "SumLargest", function(.Object, ..., x, k) {
  .Object@x <- x
  .Object@k <- k
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{SumLargest} object.
#' @param values A list of arguments to the atom.
#' @describeIn SumLargest The sum of the \code{k} largest entries of the vector or matrix.
setMethod("to_numeric", "SumLargest", function(object, values) {
  # Return the sum of the k largest entries of the matrix
  value <- as.numeric(values[[1]])
  k <- min(object@k, length(value))
  val_sort <- sort(value, decreasing = TRUE)
  sum(val_sort[1:k])
})

#' @describeIn SumLargest Check that \code{k} is a positive integer.
setMethod("validate_args",   "SumLargest", function(object) {
  if(as.integer(object@k) != object@k || object@k <= 0)
    stop("[SumLargest: validation] k must be a positive integer")
  callNextMethod()
})

#' @describeIn SumLargest The atom is a scalar.
setMethod("dim_from_args", "SumLargest", function(object) { c() })

#' @describeIn SumLargest The sign of the atom.
setMethod("sign_from_args", "SumLargest", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn SumLargest The atom is convex.
setMethod("is_atom_convex", "SumLargest", function(object) { TRUE })

#' @describeIn SumLargest The atom is not concave.
setMethod("is_atom_concave", "SumLargest", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SumLargest The atom is weakly increasing in every argument.
setMethod("is_incr", "SumLargest", function(object, idx) { TRUE })

#' @describeIn SumLargest The atom is not weakly decreasing in any argument.
setMethod("is_decr", "SumLargest", function(object, idx) { FALSE })

#' @describeIn SumLargest A list containing \code{k}.
setMethod("get_data", "SumLargest", function(object) { list(object@k) })

setMethod(".grad", "SumLargest", function(object, values) {
  # Grad: 1 for each of the k largest indices
  value <- as.numeric(t(values[[1]]))
  k <- min(object@k, length(value))
  indices <- order(value, decreasing = TRUE)
  arg_dim <- dim(object@args[[1]])
  D <- matrix(0, nrow = arg_dim[1]*arg_dim[2], ncol = 1)
  D[indices[1:k]] <- 1
  list(Matrix(D, sparse = TRUE))
})

SumSmallest <- function(x, k) {
  x <- as.Constant(x)
  -SumLargest(x = -x, k = k)
}

SumSquares <- function(expr) { QuadOverLin(x = expr, y = 1) }

TotalVariation <- function(value, ...) {
  value <- as.Constant(value)
  if(ndim(value) == 0)
    stop("TotalVariation cannot take a scalar argument")
  else if(ndim(value) == 1)   # L1 norm for vectors
    Norm(value[-1] - value[1:(nrow(value)-1)], 1)
  else {   # L2 norm for matrices
    val_dim <- dim(value)
    rows <- val_dim[1]
    cols <- val_dim[2]
    args <- lapply(list(...), as.Constant)
    values <- c(list(value), args)
    
    diffs <- list()
    for(mat in values) {
      diffs <- c(diffs, list(mat[1:(rows-1), 2:cols] - mat[1:(rows-1), 1:(cols-1)],
                             mat[2:rows, 1:(cols-1)] - mat[1:(rows-1), 1:(cols-1)]))
    }
    len <- nrow(diffs[[1]]) * ncol(diffs[[2]])
    stacked <- .VStack(atom_args = lapply(diffs, function(diff) { Reshape(diff, c(1, len)) }))
    SumEntries(Norm(stacked, p = 2, axis = 2))
  }
}
