#'
#' The Atom class.
#'
#' This virtual class represents atomic expressions in CVXR.
#'
#' @name Atom-class
#' @aliases Atom
#' @rdname Atom-class
Atom <- setClass("Atom", representation(args = "list", .size = "numeric"), prototype(args = list(), .size = NA_real_),
                 validity = function(object) {
                   if(length(object@args) == 0)
                     stop("[Atom: args] no arguments given to ", class(object))
                   return(TRUE)
                 }, contains = c("VIRTUAL", "Expression"))

setMethod("initialize", "Atom", function(.Object, ..., args = list(), .size = NA_real_) {
  .Object@args <- lapply(args, function(arg) { as.Constant(arg) })
  validate_args(.Object)
  .Object@.size <- size_from_args(.Object)
  callNextMethod(.Object, ...)
})

setMethod("show", "Atom", function(object) {
  cat(class(object), "(", paste(lapply(object@args, function(arg) { as.character(arg) }), collapse = ", "), ")", sep = "")
})

#' @param x,object An \linkS4class{Atom} object.
#' @describeIn Atom Raises an error if the arguments are invalid.
setMethod("validate_args", "Atom", function(object) { return() })

#' @rdname size_from_args
setMethod("size_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom The \code{c(row, col)} dimensions of the atom.
setMethod("size", "Atom", function(object) { object@.size })

#' @describeIn Atom The \code{c(row, col)} dimensions of the atom.
setMethod("dim", "Atom", function(x) { size(x) })

#' @describeIn Atom The number of rows in the atom.
setMethod("nrow", "Atom", function(x) { size(x)[1] })

#' @describeIn Atom The number of columns in the atom.
setMethod("ncol", "Atom", function(x) { size(x)[2] })

#' @rdname sign_from_args
setMethod("sign_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom A logical value indicating whether the atom is positive.
setMethod("is_positive", "Atom", function(object) { sign_from_args(object)[1] })

#' @describeIn Atom A logical value indicating whether the atom is negative.
setMethod("is_negative", "Atom", function(object) { sign_from_args(object)[2] })

#' @rdname curvature-atom
setMethod("is_atom_convex", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_concave", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_affine", "Atom", function(object) { is_atom_concave(object) && is_atom_convex(object) })

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

#' @describeIn Atom Represent the atom as an affine objective and conic constraints.
setMethod("canonicalize", "Atom", function(object) {
  # Constant atoms are treated as a leaf
  if(is_constant(object)) {
    # Parameterized expressions are evaluated later
    if(!is.na(parameters(object)) && length(parameters(object)) > 0) {
      size <- size(object)
      param <- CallbackParam(value(object), size[1], size[2])
      return(canonical_form(param))
    # Non-parameterized expressions are evaluated immediately
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
    # Special info required by the graph implementation
    data <- get_data(object)
    graph <- graph_implementation(object, arg_objs, size(object), data)
    return(list(graph[[1]], c(constraints, graph[[2]])))
  }
})

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Atom The graph implementation of the atom.
setMethod("graph_implementation", "Atom", function(object, arg_objs, size, data = NA_real_) { stop("Unimplemented") })

#' @describeIn Atom List of \linkS4class{Variable} objects in the atom.
setMethod("variables", "Atom", function(object) {
  var_list <- lapply(object@args, function(arg) { variables(arg) })
  unique(flatten_list(var_list))
})

#' @describeIn Atom List of \linkS4class{Parameter} objects in the atom.
setMethod("parameters", "Atom", function(object) {
  param_list <- lapply(object@args, function(arg) { parameters(arg) })
  unique(flatten_list(param_list))
})

#' @describeIn Atom List of \linkS4class{Constant} objects in the atom.
setMethod("constants", "Atom", function(object) {
  const_list <- lapply(object@args, function(arg) { constants(arg) })
  unique(flatten_list(const_list))   # TODO: Is this the correct way to remove duplicates?
})

#' @describeIn Atom The value of the atom.
setMethod("value", "Atom", function(object) {
  # Catch the case when the expression is known to be zero through DCP analysis
  if(is_zero(object)) {
    size <- size(object)
    result <- matrix(0, nrow = size[1], ncol = size[2])
  } else {
    arg_values <- list()
    idx <- 1
    for(arg in object@args) {
      # An argument without a value makes all higher level values NA.
      # But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
      arg_val <- value(arg)
      if(any(is.na(arg_val)) && !is_constant(object))
        return(NA)
      else {
        arg_values[[idx]] <- arg_val
        idx <- idx + 1
      }
    }
    result <- to_numeric(object, arg_values)
  }

  # Reduce to scalar if possible
  if(all(intf_size(result) == c(1,1)))
    intf_scalar_value(result)
  else
    result
})

setMethod(".grad", "Atom", function(object, values) { stop("Unimplemented") })

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
    # A dictionary of gradients wrt variables
    # Partial argument / partial x
    grad_arg <- grad(arg)
    for(key in names(grad_arg)) {
      # None indicates gradient is not defined
      if(any(is.na( as.numeric(grad_arg[[key]]) )) || any(is.na( as.numeric(grad_self[[idx]]) )))
        result[[key]] <- NA_real_
      else {
        D <- grad_arg[[key]] %*% grad_self[[idx]]
        # Convert 1x1 matrices to scalars
        if((is.matrix(D) || is(D, "Matrix")) && dim(D) == c(1,1))
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

setMethod(".domain", "Atom", function(object) { list() })

#' @describeIn Atom A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "Atom", function(object) {
  cons <- list()
  for(arg in object@args) {
    for(con in domain(arg))
      cons <- c(cons, con)
  }
  c(.domain(object), cons)
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
AxisAtom <- setClass("AxisAtom", representation(expr = "ConstValORExpr", axis = "ANY"), prototype(axis = NA_real_), 
                     validity = function(object) {
                       if(length(object@axis) != 1 || !(is.numeric(axis) || is.logical(axis)))
                         stop("[AxisAtom: axis] axis must equal 1 (row), 2 (column), or NA (row and column)")
                       return(TRUE)
                     }, contains = c("VIRTUAL", "Atom"))

setMethod("initialize", "AxisAtom", function(.Object, ..., expr, axis) {
  .Object@expr <- expr
  .Object@axis <- axis
  .Object <- callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @param object An \linkS4class{Atom} object.
#' @describeIn AxisAtom The size of the atom deteremined from its arguments.
setMethod("size_from_args", "AxisAtom", function(object) {
  if(is.na(object@axis))
    c(1, 1)
  else if(object@axis == 1)
    c(size(object@args[[1]])[1], 1)
  else   # axis == 2
    c(1, size(object@args[[1]])[2])
})

#' @describeIn AxisAtom A list containing \code{axis}.
setMethod("get_data", "AxisAtom", function(object) { list(object@axis) })

#' @describeIn AxisAtom Check that the new shape has the same number of entries as the old.
setMethod("validate_args", "AxisAtom", function(object) {
  if(length(object@axis) != 1 || !(is.na(object@axis) || object@axis %in% c(1,2)))
     stop("Invalid argument for axis: must equal 1 (row), 2 (column), or NA (row and column)")
})

setMethod(".axis_grad", "AxisAtom", function(object, values) {
  m <- size(object@args[[1]])[1]
  n <- size(object@args[[1]])[2]
  if(is.na(object@axis)) {
    value <- matrix(values[[1]], nrow = m*n, ncol = 1)
    D <- .column_grad(object, value)
    if(is(D, "Matrix") || !any(is.na(D)))
      D <- Matrix(D, sparse = TRUE)
  } else {
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
#' The AffineProd class.
#'
#' This class represents the product of two affine expressions.
#'
#' @slot x An \linkS4class{Expression} or numeric constant representing the left-hand value.
#' @slot y An \linkS4class{Expression} or numeric constant representing the right-hand value.
#' @name AffineProd-class
#' @aliases AffineProd
#' @rdname AffineProd-class
.AffineProd <- setClass("AffineProd", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric constant representing the left-hand value.
#' @param y An \linkS4class{Expression} or numeric constant representing the right-hand value.
#' @rdname AffineProd-class
AffineProd <- function(x, y) { .AffineProd(x = x, y = y) }

setMethod("initialize", "AffineProd", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., args = list(x, y))
})

#' @describeIn AffineProd Check dimensions of arguments and linearity.
setMethod("validate_args", "AffineProd", function(object) {
  if(!is_affine(object@args[[1]]) || !is_affine(object@args[[2]]))
    stop("The arguments to AffineProd must be affine")
  mul_shapes(size(object@args[[1]]), size(object@args[[2]]))
})

#' @param object An \linkS4class{AffineProd} object.
#' @param values A list of arguments to the atom.
#' @describeIn AffineProd The product of two affine expressions.
setMethod("to_numeric", "AffineProd", function(object, values) { values[[1]] %*% values[[2]] })

#' @describeIn AffineProd The size of the atom.
setMethod("size_from_args", "AffineProd", function(object) { mul_shapes(size(object@args[[1]]), size(object@args[[2]])) })

#' @describeIn AffineProd Default to rules for times.
setMethod("sign_from_args", "AffineProd", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @describeIn AffineProd Affine times affine is not convex.
setMethod("is_atom_convex", "AffineProd", function(object) { FALSE })

#' @describeIn AffineProd Affine times affine is not concave.
setMethod("is_atom_concave", "AffineProd", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn AffineProd A logical value indicating whether the atom is weakly increasing in \code{idx}.
setMethod("is_incr", "AffineProd", function(object, idx) { is_positive(object@args[[2-idx]]) })

#' @describeIn AffineProd A logical value indicating whether the atom is weakly decreasing in \code{idx}.
setMethod("is_decr", "AffineProd", function(object, idx) { is_negative(object@args[[2-idx]]) })

#' @describeIn AffineProd The affine product is always quadratic.
setMethod("is_quadratic", "AffineProd", function(object) { TRUE })

setMethod(".grad", "AffineProd", function(object, values) {
  X <- values[[1]]
  Y <- values[[2]]
  size11 <- size(object@args[[1]])[1]
  size22 <- size(object@args[[2]])[2]

  DX_rows <- prod(size(object@args[[1]]))
  cols <- size11 * size22

  # DX = [diag(Y11), diag(Y12), ...]
  #      [diag(Y21), diag(Y22), ...]
  #      [  ...        ...      ...]
  if(!(is.matrix(Y) || is(Y, "Matrix")))
    Y <- as.matrix(Y)
  DX <- do.call("rbind", apply(Y, 1, function(row) {
      do.call("cbind", lapply(row, function(x) { sparseMatrix(i = 1:size11, j = 1:size11, x = x) }))
    }))
  DY <- bdiag(lapply(1:size22, function(i) { t(X) }))
  list(DX, DY)
})

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
  .Object <- callNextMethod(.Object, ..., args = list(.Object@x))

  x <- .Object@args[[1]]
  if(size(x)[1] == 1)
    n <- size(x)[2]
  else if(size(x)[2] == 1)
    n <- size(x)[1]
  else
    stop("x must be a row or column vector")

  if(any(is.na(p)))
    p <- rep(1, n)

  if(length(p) != n)
    stop("x and p must have the same number of elements")

  if(any(p < 0) || sum(p) <= 0)
    stop("powers must be nonnegative and not all zero")

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

#' @describeIn GeoMean Empty function since validation of arguments is done during atom initialization.
setMethod("validate_args", "GeoMean", function(object) { return() })

#' @param object A \linkS4class{GeoMean} object.
#' @param values A list of arguments to the atom.
#' @describeIn GeoMean The (weighted) geometric mean of the elements of \code{x}.
setMethod("to_numeric", "GeoMean", function(object, values) {
  values <- as.numeric(values[[1]])

  if(requireNamespace("Rmpfr") && requireNamespace("gmp")) {
    val <- 1.0
    for(idx in 1:length(values)) {
      x <- values[[idx]]
      p <- object@w[idx]
      val <- val * Rmpfr::mpfr(x, Rmpfr::getPrec(x))^p
    }
    return(gmp::asNumeric(val))   # TODO: Handle mpfr objects in the backend later
  }

  val <- mapply(function(x, p) { x^p }, values, gmp::asNumeric(object@w))
  Reduce("*", val)
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
setMethod("size_from_args", "GeoMean", function(object) { c(1,1) })

#' @describeIn GeoMean The atom is non-negative.
setMethod("sign_from_args", "GeoMean", function(object) { c(TRUE, FALSE) })

#' @describeIn GeoMean The atom is not convex.
setMethod("is_atom_convex", "GeoMean", function(object) { FALSE })

#' @describeIn GeoMean The atom is concave.
setMethod("is_atom_concave", "GeoMean", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn GeoMean The atom is weakly increasing in every argument.
setMethod("is_incr", "GeoMean", function(object, idx) { TRUE })

#' @describeIn GeoMean The atom is not weakly decreasing in any argument.
setMethod("is_decr", "GeoMean", function(object, idx) { FALSE })

#' @describeIn GeoMean Returns \code{list(w, dyadic completion, tree of dyads)}.
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

GeoMean.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  w <- data[[1]]
  w_dyad <- data[[2]]
  tree <- data[[3]]
  t <- create_var(c(1,1))

  if(size(arg_objs[[1]])[2] == 1)
    x_list <- lapply(1:length(w), function(i) { Index.get_index(arg_objs[[1]], list(), i, 1)$idx })
  else if(size(arg_objs[[1]])[1] == 1)
    x_list <- lapply(1:length(w), function(i) { Index.get_index(arg_objs[[1]], list(), 1, i)$idx })

  # TODO: Catch cases where we have (0,0,1)?
  # TODO: What about curvature (should be affine) in trivial case of (0,0,1),
  # should this behavior match what we do in power?
  list(t, gm_constrs(t, x_list, w))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn GeoMean The graph implementation of the atom.
setMethod("graph_implementation", "GeoMean", function(object, arg_objs, size, data = NA_real_) {
  GeoMean.graph_implementation(arg_objs, size, data)
})

HarmonicMean <- function(x) {
  x <- as.Constant(x)
  prod(size(x)) * Pnorm(x = x, p = -1)
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
  callNextMethod(.Object, ..., args = list(.Object@A))
})

#' @describeIn LambdaMax Check that \code{A} is square.
setMethod("validate_args", "LambdaMax", function(object) {
  if(size(object@args[[1]])[1] != size(object@args[[1]])[2])
    stop("The argument to LambdaMax must resolve to a square matrix")
})

#' @param object A \linkS4class{LambdaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn LambdaMax The largest eigenvalue of \code{A}. Requires that \code{A} be symmetric.
setMethod("to_numeric", "LambdaMax", function(object, values) {
  if(!all(t(values[[1]]) == values[[1]]))
    stop("LambdaMax called on a non-symmetric matrix")
  max(eigen(values[[1]], only.values = TRUE)$values)
})

#' @describeIn LambdaMax The atom is a scalar.
setMethod("size_from_args", "LambdaMax", function(object) { c(1, 1) })

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

setMethod(".domain", "LambdaMax", function(object) { list(t(object@args[[1]]) == object@args[[1]]) })

setMethod(".grad", "LambdaMax", function(object, values) {
  r <- eigen(values[[1]], only.values = FALSE)
  v <- r$vectors  # eigenvectors
  w <- r$values   # eigenvalues

  d <- rep(0, length(w))
  d[1] <- 1
  d <- diag(d)
  D <- v %*% d %*% t(v)
  list(Matrix(as.numeric(D), sparse = TRUE))
})

LambdaMax.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]
  n <- size(A)[1]
  # SDP constraint.
  t <- create_var(c(1,1))
  prom_t <- lo.promote(t, c(n,1))
  # I*t - A
  expr <- lo.sub_expr(lo.diag_vec(prom_t), A)
  list(t, list(SDP(expr)))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn LambdaMax The graph implementation of the atom.
setMethod("graph_implementation", "LambdaMax", function(object, arg_objs, size, data = NA_real_) {
  LambdaMax.graph_implementation(arg_objs, size, data)
})

LambdaMin <- function(A) {
  A <- as.Constant(A)
  -LambdaMax(-A)
}

LambdaSumLargest <- function(A, k) {
  A <- as.Constant(A)
  if(size(A)[1] != size(A)[2])
    stop("First argument must be a square matrix")
  else if(as.integer(k) != k || k <= 0)
    stop("Second argument must be a positive integer")

  Z <- Semidef(size(A)[1])
  k*LambdaMax(A - Z) + Trace(Z)
}

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
  callNextMethod(.Object, ..., args = list(.Object@A))
})

#' @describeIn LogDet Check that \code{A} is square.
setMethod("validate_args", "LogDet", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("The argument to LogDet must be a square matrix")
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

#' @describeIn LogDet The atom is a scalar.
setMethod("size_from_args", "LogDet", function(object) { c(1, 1) })

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

LogDet.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]  # n by n matrix
  n <- size(A)[1]
  X <- create_var(c(2*n, 2*n))
  canon <- canonical_form(Semidef(2*n))
  X <- canon[[1]]
  constraints <- canon[[2]]
  Z <- create_var(c(n,n))
  D <- create_var(c(n,1))

  # Require that X and A are PSD
  constraints <- c(constraints, list(SDP(A)))

  # Fix Z as upper triangular, D as diagonal, and diag(D) as diag(Z)
  Z_lower_tri <- lo.upper_tri(lo.transpose(Z))
  constraints <- c(constraints, list(create_eq(Z_lower_tri)))

  # D[i,i] = Z[i,i]
  constraints <- c(constraints, list(create_eq(D, lo.diag_mat(Z))))

  # Fix X using the fact that A must be affine by the DCP rules
  # X[1:n, 1:n] == D
  constraints <- Index.block_eq(X, lo.diag_vec(D), constraints, 1, n, 1, n)

  # X[1:n, (n+1):(2*n)] == Z
  constraints <- Index.block_eq(X, Z, constraints, 1, n, n+1, 2*n)

  # X[(n+1):(2*n), (n+1):(2*n)] == A
  constraints <- Index.block_eq(X, A, constraints, n+1, 2*n, n+1, 2*n)

  # Add the objective sum(log(D[i,i]))
  graph <- Log.graph_implementation(list(D), c(n, 1))
  obj <- graph[[1]]
  constr <- graph[[2]]
  list(lo.sum_entries(obj), c(constraints, constr))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn LogDet The graph implementation of the atom.
setMethod("graph_implementation", "LogDet", function(object, arg_objs, size, data = NA_real_) {
  LogDet.graph_implementation(arg_objs, size, data)
})

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

#' @describeIn LogSumExp The atom is positive.
setMethod("sign_from_args",  "LogSumExp", function(object) { c(TRUE, FALSE) })

#' @describeIn LogSumExp The atom is convex.
setMethod("is_atom_convex", "LogSumExp", function(object) { TRUE })

#' @describeIn LogSumExp The atom is not concave.
setMethod("is_atom_concave", "LogSumExp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is not monotonic in any argument.
setMethod("is_incr", "LogSumExp", function(object, idx) { FALSE })

#' @describeIn LogSumExp The atom is not monotonic in any argument.
setMethod("is_decr", "LogSumExp", function(object, idx) { FALSE })

LogSumExp.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  axis <- data[[1]]
  t <- create_var(size)

  # sum(exp(x - t)) <= 1
  if(is.na(axis)) {
    prom_t <- lo.promote(t, size(x))
    expr <- lo.sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]
    obj <- lo.sum_entries(obj)
  } else if(axis == 2) {
    prom_size <- c(size(x)[1], 1)
    ones <- create_const(matrix(1, nrow = prom_size[1], ncol = prom_size[2]), prom_size)
    prom_t <- lo.mul_expr(ones, t, size(x))
    expr <- lo.sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]

    const_size <- c(1, size(x)[1])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- lo.mul_expr(ones, obj, size)
  } else {    # axis == 1
    prom_size <- c(1, size(x)[2])
    ones <- create_const(matrix(1, nrow = prom_size[1], ncol = prom_size[2]), prom_size)
    prom_t <- lo.rmul_expr(t, ones, size(x))
    expr <- lo.sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]

    const_size <- c(size(x)[2], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- lo.rmul_expr(obj, ones, size)
  }

  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  constraints <- c(constraints, list(create_leq(obj, ones)))
  list(t, constraints)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn LogSumExp The graph implementation of the atom.
setMethod("graph_implementation", "LogSumExp", function(object, arg_objs, size, data = NA_real_) {
  LogSumExp.graph_implementation(arg_objs, size, data)
})

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
MatrixFrac <- function(X, P) { .MatrixFrac(X = X, P = P) }

setMethod("initialize", "MatrixFrac", function(.Object, ..., X, P) {
  .Object@X <- X
  .Object@P <- P
  callNextMethod(.Object, ..., args = list(.Object@X, .Object@P))
})

#' @describeIn MatrixFrac Check that the dimensions of \code{x} and \code{P} match.
setMethod("validate_args", "MatrixFrac", function(object) {
  X <- object@args[[1]]
  P <- object@args[[2]]
  if(size(P)[1] != size(P)[2])
    stop("The second argument to MatrixFrac must be a square matrix")
  else if(size(X)[1] != size(P)[1])
    stop("The arguments to MatrixFrac have incompatible dimensions")
})

#' @param object A \linkS4class{MatrixFrac} object.
#' @param values A list of arguments to the atom.
#' @describeIn MatrixFrac The trace of \eqn{X^TP^{-1}X}.
setMethod("to_numeric", "MatrixFrac", function(object, values) {
  # TODO: Raise error if not invertible?
  X <- values[[1]]
  P <- values[[2]]
  sum(diag(t(X) %*% base::solve(P) %*% X))
})

#' @describeIn MatrixFrac The atom is a scalar.
setMethod("size_from_args", "MatrixFrac", function(object) { c(1, 1) })

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

setMethod(".domain", "MatrixFrac", function(object) { list(object@args[[2]] %>>% 0) })

setMethod(".grad", "MatrixFrac", function(object, values) {
  X <- as.matrix(values[[1]])
  P <- as.matrix(values[[2]])
  P_inv <- tryCatch({
    base::solve(P)
  }, error = function(e) {
    return(NA_real_)
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

MatrixFrac.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  X <- arg_objs[[1]]   # n by m matrix
  P <- arg_objs[[2]]   # n by n matrix
  size <- size(X)
  n <- size[1]
  m <- size[2]

  # Create a matrix with Schur complement Tmat - t(X) * P^-1 * X
  M <- create_var(c(n+m, n+m))
  Tmat <- create_var(c(m, m))
  constraints <- list()

  # Fix M using the fact that P must be affine by the DCP rules.
  # M[1:n, 1:n] == P
  constraints <- Index.block_eq(M, P, constraints, 1, n, 1, n)

  # M[1:n, (n+1):(n+m)] == X
  constraints <- Index.block_eq(M, X, constraints, 1, n, n+1, n+m)

  # M[(n+1):(n+m), (n+1):(n+m)] == Tmat
  constraints <- Index.block_eq(M, Tmat, constraints, n+1, n+m, n+1, n+m)

  # Add SDP constraints.
  list(lo.trace(Tmat), c(constraints, list(SDP(M))))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn MatrixFrac The graph implementation of the atom.
setMethod("graph_implementation", "MatrixFrac", function(object, arg_objs, size, data = NA_real_) {
  MatrixFrac.graph_implementation(arg_objs, size, data)
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
MaxEntries <- function(x, axis = NA_real_) { .MaxEntries(expr = x, axis = axis) }

#' @param object A \linkS4class{MaxEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MaxEntries The largest entry in \code{x}.
setMethod("to_numeric", "MaxEntries", function(object, values) {
  if(is.na(object@axis))
    max(values[[1]])
  else
    apply(values[[1]], object@axis, max)
})

#' @describeIn MaxEntries The sign of the atom.
setMethod("sign_from_args",  "MaxEntries", function(object) { c(is_positive(object@args[[1]]), is_negative(object@args[[1]])) })

#' @describeIn MaxEntries The atom is convex.
setMethod("is_atom_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries The atom is not concave.
setMethod("is_atom_concave", "MaxEntries", function(object) { FALSE })

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
  idx <- which.max(value)
  D <- rep(0, length(value))
  D[idx] <- 1
  D
})

MaxEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  axis <- data[[1]]
  if(is.na(axis)) {
    t <- create_var(c(1,1))
    promoted_t <- lo.promote(t, size(arg_objs[[1]]))
  } else if(axis == 2) {
    t <- create_var(c(1, size(arg_objs[[1]])[2]))
    const_size <- c(size(arg_objs[[1]])[1], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    promoted_t <- lo.mul_expr(ones, t, size(arg_objs[[1]]))
  } else {   # axis == 1
    t <- create_var(c(size(arg_objs[[1]])[1], 1))
    const_size <- c(1, size(arg_objs[[1]])[2])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    promoted_t <- lo.rmul_expr(t, ones, size(arg_objs[[1]]))
  }

  constraints <- list(create_leq(arg_objs[[1]], promoted_t))
  list(t, constraints)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn MaxEntries The graph implementation of the atom.
setMethod("graph_implementation", "MaxEntries", function(object, arg_objs, size, data = NA_real_) {
  MaxEntries.graph_implementation(arg_objs, size, data)
})

MinEntries <- function(x, axis = NA_real_) {
  x <- as.Constant(x)
  -MaxEntries(-x, axis = axis)
}

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
.Pnorm <- setClass("Pnorm", representation(p = "numeric", max_denom = "numeric", .approx_error = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_), contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param p A number greater than or equal to 1, or equal to positive infinity.
#' @param max_denom The maximum denominator considered in forming a rational approximation for \eqn{p}.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @rdname Pnorm-class
Pnorm <- function(x, p = 2, axis = NA_real_, max_denom = 1024) { .Pnorm(expr = x, axis = axis, p = p, max_denom = max_denom) }

setMethod("initialize", "Pnorm", function(.Object, ..., p = 2, max_denom = 1024, .approx_error = NA_real_) {
  .Object@max_denom <- max_denom

  # TODO: Deal with fractional powers correctly
  p_old <- p
  # if(p == Inf)
  #  .Object@p <- Inf
  # else if(p < 0)
  #  .Object@p <- pow_neg(p, max_denom)
  # else if(p > 0 && p < 1)
  #  .Object@p <- pow_mid(p, max_denom)
  # else if(p > 1)
  #  .Object@p <- pow_high(p, max_denom)
  # else if(p == 1)
  #  .Object@p <- 1
  # else
  #  stop("[Pnorm: validation] Invalid value for p ", p)

  .Object@p <- p
  if(.Object@p == Inf)
    .Object@.approx_error <- 0
  else
    .Object@.approx_error <- abs(.Object@p - p_old)
  callNextMethod(.Object, ...)
})

#' @describeIn Pnorm Check that the arguments are valid.
setMethod("validate_args", "Pnorm", function(object) {
  callNextMethod()
  if(!is.na(object@axis) && object@p != 2)
    stop("The axis parameter is only supported for p = 2")
})

#' @describeIn Pnorm The name and arguments of the atom.
#' @export
setMethod("name", "Pnorm", function(object) {
  sprintf("%s(%s, %s)", class(object), name(object@args[1]), object@p)
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

  if(is.na(object@axis))
    retval <- .p_norm(values, object@p)
  else
    retval <- apply(values, object@axis, function(x) { .p_norm(x, object@p) })
  retval
})

#' @describeIn Pnorm The atom is positive.
setMethod("sign_from_args",  "Pnorm", function(object) { c(TRUE, FALSE) })

#' @describeIn Pnorm The atom is convex if \eqn{p \geq 1}.
setMethod("is_atom_convex", "Pnorm", function(object) { object@p >= 1})

#' @describeIn Pnorm The atom is concave if \eqn{p < 1}.
setMethod("is_atom_concave", "Pnorm", function(object) { object@p < 1 })

#' @param idx An index into the atom.
#' @describeIn Pnorm The atom is weakly increasing if \eqn{p < 1} or \eqn{p \geq 1} and \code{x} is positive.
setMethod("is_incr", "Pnorm", function(object, idx) { object@p < 1 || (object@p >= 1 && is_positive(object@args[[1]])) })

#' @describeIn Pnorm The atom is weakly decreasing if \eqn{p \geq 1} and \code{x} is negative.
setMethod("is_decr", "Pnorm", function(object, idx) { object@p >= 1 && is_negative(object@args[[1]]) })

#' @describeIn Pnorm The atom is piecewise linear only if \code{x} is piecewise linear, and either \eqn{p = 1} or \eqn{p = \infty}.
setMethod("is_pwl", "Pnorm", function(object) { (object@p == 1 || object@p == Inf) && is_pwl(object@args[[1]]) })

#' @describeIn Pnorm Returns \code{list(p, axis)}.
setMethod("get_data", "Pnorm", function(object) { list(object@p, object@axis) })

setMethod(".domain", "Pnorm", function(object) {
  if(object@p < 1 && object@p != 0)
    list(object@args[[1]] >= 0)
  else
    list()
})

setMethod(".grad", "Pnorm", function(object, values) { .axis_grad(object, values) })

setMethod(".column_grad", "Pnorm", function(object, value) {
  rows <- prod(size(object@args[[1]]))
  value <- as.matrix(value)

  # Outside domain
  if(object@p < 1 && any(value <= 0))
    return(NA_real_)
  D_null <- sparseMatrix(i = c(), j = c(), dims = c(rows, 1))
  if(object@p == 1) {
    D_null <- D_null + (value > 0)
    D_null <- D_null - (value < 0)
    return(Matrix(as.numeric(D_null), sparse = TRUE))   # TODO: Is this redundant? Check against CVXPY
  }
  denominator <- .p_norm(value, object@p)
  denominator <- denominator^(object@p - 1)

  # Subgrad is 0 when denom is 0 (or undefined)
  if(denominator == 0) {
    if(object@p >= 1)
      return(D_null)
    else
      return(NA_real_)
  } else {
    numerator <- value^(object@p - 1)
    frac <- numerator / denominator
    return(matrix(as.numeric(frac)))
  }
})

Pnorm.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  p <- data[[1]]
  axis <- data[[2]]
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  constraints <- list()

  # First, take care of special cases p = 2, Inf, and 1
  if(p == 2) {
    if(is.na(axis))
      return(list(t, list(SOC(t, list(x)))))
    else {
      t <- create_var(size)
      return(list(t, list(SOCAxis(lo.reshape(t, c(prod(t$size), 1)), x, axis))))
    }
  }

  if(p == Inf) {
    t_ <- lo.promote(t, x$size)
    return(list(t, list(create_leq(x, t_), create_geq(lo.sum_expr(list(x, t_))))))
  }

  # We need absolute value constraint for symmetric convex branches (p >= 1)
  # We alias |x| as x from this point forward to make code pretty
  if(p >= 1) {
    absx <- create_var(x$size)
    constraints <- c(constraints, list(create_leq(x, absx), create_geq(lo.sum_expr(list(x, absx))) ))
    x <- absx
  }

  if(p == 1)
    return(list(lo.sum_entries(x), constraints))

  # Now take care of remaining convex/concave branches
  # To create rational powers, need new variable r and constraint sum(r) == t
  r <- create_var(x$size)
  t_ <- lo.promote(t, x$size)
  constraints <- c(constraints, list(create_eq(lo.sum_entries(r), t)))

  p <- as.bigq(p)   # TODO: Can we simplify the fraction, e.g. for p = 1.6?
  if(p < 0)
    constraints <- c(constraints, gm_constrs(t_, list(x, r), c(-p/(1-p), 1/(1-p)) ))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, t_), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, t_), c(1/p, 1-1/p)))

  list(t, constraints)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Pnorm The graph implementation of the atom.
setMethod("graph_implementation", "Pnorm", function(object, arg_objs, size, data = NA_real_) {
  Pnorm.graph_implementation(arg_objs, size, data)
})

Norm <- function(x, p = 2, axis = NA_real_) {
  x <- as.Constant(x)

  # Norms for scalars same as absolute value
  if(p == 1 || is_scalar(x))
    Pnorm(x = x, p = 1, axis = axis)
  else if(p == "inf" || p == Inf)
    Pnorm(x = x, p = Inf, axis = axis)
  else if(p == "nuc")
    NormNuc(A = x)
  else if(p == "fro")
    Pnorm(x = x, p = 2, axis = axis)
  else if(p == 2) {
    if(is.na(axis) && is_matrix(x))
      SigmaMax(A = x)
    else
      Pnorm(x = x, p = 2, axis = axis)
  } else
    Pnorm(x = x, p = p, axis = axis)
}

Norm1   <- function(x, axis = NA_real_) { Pnorm(x = x, p = 1, axis = axis) }
Norm2   <- function(x, axis = NA_real_) { Pnorm(x = x, p = 2, axis = axis) }
NormInf <- function(x, axis = NA_real_) { Pnorm(x = x, p = Inf, axis = axis) }
MixedNorm <- function(X, p = 2, q = 1) {
  X <- as.Constant(X)

  # Inner norms
  vecnorms <- lapply(1:size(X)[1], function(i) { Norm(X[i,], p) })

  # Outer norms
  Norm(new("HStack", args = vecnorms), q)
}

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
  callNextMethod(.Object, ..., args = list(.Object@A))
})

#' @param object A \linkS4class{NormNuc} object.
#' @param values A list of arguments to the atom.
#' @describeIn NormNuc The nuclear norm (i.e., the sum of the singular values) of \code{A}.
setMethod("to_numeric", "NormNuc", function(object, values) {
  # Returns the nuclear norm (i.e. the sum of the singular values) of A
  sum(svd(values[[1]])$d)
})

#' @describeIn NormNuc The atom is a scalar.
setMethod("size_from_args", "NormNuc", function(object) { c(1, 1) })

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
  list(Matrix(as.numeric(D), sparse = TRUE))   # TODO: Make sure D is vectorized correctly
})

NormNuc.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]
  size <- size(A)
  rows <- size[1]
  cols <- size[2]

  # Create the equivalent problem:
  # minimize (trace(U) + trace(V))/2
  # subject to: [U A; t(A) V] is positive semidefinite
  X <- create_var(c(rows+cols, rows+cols))
  constraints <- list()

  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:rows, (rows+1):(rows+cols)] == A
  constraints <- Index.block_eq(X, A, constraints, 1, rows, rows+1, rows+cols)
  half <- create_const(0.5, c(1,1))
  trace_expr <- lo.mul_expr(half, lo.trace(X), c(1, 1))

  # Add SDP constraint.
  list(trace_expr, c(list(SDP(X)), constraints))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn NormNuc The graph implementation of the atom.
setMethod("graph_implementation", "NormNuc", function(object, arg_objs, size, data = NA_real_) {
  NormNuc.graph_implementation(arg_objs, size, data)
})

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

  # Check dimensions
  n <- size(P)[1]
  if(size(P)[2] != n || any(size(x) != c(n,1)))
    stop("Invalid dimensions for arguments")
  if(length(parameters(P)) > 0)
    stop("P cannot be a parameter")
  if(is_constant(x))
    return(t(x) %*% P %*% x)
  else if(is_constant(P)) {
    P <- as.matrix(value(P))

    # Force symmetry
    P <- (P + t(P)) / 2.0
    decomp <- .decomp_quad(P)
    scale <- decomp[[1]]
    M1 <- decomp[[2]]
    M2 <- decomp[[3]]

    ret <- 0
    if(prod(dim(M1)) > 0)
      ret <- ret + scale * SumSquares(Constant(t(M1)) %*% x)
    if(prod(dim(M2)) > 0)
      ret <- ret - scale * SumSquares(Constant(t(M2)) %*% x)
    return(ret)
  } else
    stop("At least one argument to QuadForm must be constant")
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
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y))
})

#' @describeIn QuadOverLin Check the dimensions of the arguments.
setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@args[[2]]))
    stop("[QuadOverLin: validation] y must be a scalar")
})

#' @param object A \linkS4class{QuadOverLin} object.
#' @param values A list of arguments to the atom.
#' @describeIn QuadOverLin The sum of the entries of \code{x} squared over \code{y}.
setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(values[[1]]^2) / values[[2]] })

#' @describeIn QuadOverLin The atom is a scalar.
setMethod("size_from_args", "QuadOverLin", function(object) { c(1, 1) })

#' @describeIn QuadOverLin The atom is positive.
setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })

#' @describeIn QuadOverLin The atom is convex.
setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin The atom is not concave.
setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly increasing.
setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_positive(object@args[[idx]]) })

#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly decreasing.
setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_negative(object@args[[idx]])) || (idx == 2) })

#' @describeIn QuadOverLin True if \code{x} is affine and \code{y} is constant.
setMethod("is_quadratic", "QuadOverLin", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

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

QuadOverLin.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  y <- arg_objs[[2]]   # Known to be a scalar.
  v <- create_var(c(1,1))
  two <- create_const(2, c(1,1))
  constraints <- list(SOC(lo.sum_expr(list(y, v)),
                          list(lo.sub_expr(y, v),
                               lo.mul_expr(two, x, x$size))),
                      create_geq(y))
  list(v, constraints)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn QuadOverLin The graph implementation of the atom.
setMethod("graph_implementation", "QuadOverLin", function(object, arg_objs, size, data = NA_real_) {
  QuadOverLin.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@A))
})

#' @param object A \linkS4class{SigmaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn SigmaMax The largest singular value of \code{A}.
setMethod("to_numeric", "SigmaMax", function(object, values) { base::norm(values[[1]], type = "2") })

#' @describeIn SigmaMax The atom is a scalar.
setMethod("size_from_args", "SigmaMax", function(object) { c(1, 1) })

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

SigmaMax.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]   # n by m matrix
  size <- A$size
  n <- size[1]
  m <- size[2]

  # Create a matrix with Schur complement I*t - (1/t)*t(A)*A
  X <- create_var(c(n+m, n+m))
  t <- create_var(c(1,1))
  constraints <- list()

  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:n, 1:n] == I_n*t
  prom_t <- lo.promote(t, c(n,1))
  constraints <- Index.block_eq(X, lo.diag_vec(prom_t), constraints, 1, n, 1, n)

  # X[1:n, (n+1):(n+m)] == A
  constraints <- Index.block_eq(X, A, constraints, 1, n, n+1, n+m)

  # X[(n+1):(n+m), (n+1):(n+m)] == I_m*t
  prom_t <- lo.promote(t, c(m,1))
  constraints <- Index.block_eq(X, lo.diag_vec(prom_t), constraints, n+1, n+m, n+1, n+m)

  # Add SDP constraint.
  list(t, c(constraints, list(SDP(X))))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn SigmaMax The graph implementation of the atom.
setMethod("graph_implementation", "SigmaMax", function(object, arg_objs, size, data = NA_real_) {
  SigmaMax.graph_implementation(arg_objs, size, data)
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
.SumLargest <- setClass("SumLargest", representation(x = "ConstValORExpr", k = "numeric"),
                       validity = function(object) {
                         if(as.integer(object@k) != object@k || object@k <= 0)
                           stop("[SumLargest: validation] k must be a positive integer")
                         return(TRUE)
                         }, contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param k The number of largest values to sum over.
#' @rdname SumLargest-class
SumLargest <- function(x, k) { .SumLargest(x = x, k = k) }

setMethod("initialize", "SumLargest", function(.Object, ..., x, k) {
  .Object@x <- x
  .Object@k <- k
  callNextMethod(.Object, ..., args = list(.Object@x))
})

#' @describeIn SumLargest Check that \code{k} is a positive integer.
setMethod("validate_args",   "SumLargest", function(object) {
  if(as.integer(object@k) != object@k || object@k <= 0)
    stop("[SumLargest: validation] k must be a positive integer")
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

#' @describeIn SumLargest The atom is a scalar.
setMethod("size_from_args", "SumLargest", function(object) { c(1, 1) })

#' @describeIn SumLargest The sign of the atom.
setMethod("sign_from_args", "SumLargest", function(object) { c(is_positive(object@args[[1]]), is_negative(object@args[[1]])) })

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
  D <- matrix(0, nrow = prod(size(object@args[[1]])), ncol = 1)
  D[indices[1:k]] <- 1
  list(Matrix(D, sparse = TRUE))
})

SumLargest.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # min SumEntries(t) + k*q
  # s.t. x <= t + q, t >= 0
  x <- arg_objs[[1]]
  k <- create_const(data[[1]], c(1,1))
  q <- create_var(c(1,1))
  t <- create_var(x$size)

  sum_t <- lo.sum_entries(t)
  obj <- lo.sum_expr(list(sum_t, lo.mul_expr(k, q, c(1,1))))
  prom_q <- lo.promote(q, x$size)

  constr <- list(create_leq(x, lo.sum_expr(list(t, prom_q))), create_geq(t))
  list(obj, constr)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn SumLargest The graph implementation of the atom.
setMethod("graph_implementation", "SumLargest", function(object, arg_objs, size, data = NA_real_) {
  SumLargest.graph_implementation(arg_objs, size, data)
})

SumSmallest <- function(x, k) {
  x <- as.Constant(x)
  -SumLargest(x = -x, k = k)
}

SumSquares <- function(expr) { QuadOverLin(x = expr, y = 1) }

TotalVariation <- function(value, ...) {
  value <- as.Constant(value)
  val_size <- size(value)
  rows <- val_size[1]
  cols <- val_size[2]

  if(is_scalar(value))
    stop("TotalVariation cannot take a scalar argument")
  else if(is_vector(value))   # L1 norm for vectors
    Pnorm(value[-1] - value[1:(max(rows, cols)-1)], 1)
  else {   # L2 norm for matrices
    args <- lapply(list(...), function(arg) { as.Constant(arg) })
    values <- c(list(value), args)
    diffs <- list()
    for(mat in values) {
      diffs <- c(diffs, list(mat[1:(rows-1), 2:cols] - mat[1:(rows-1), 1:(cols-1)],
                             mat[2:rows, 1:(cols-1)] - mat[1:(rows-1), 1:(cols-1)]))
    }
    length <- size(diffs[[1]])[1] * size(diffs[[2]])[2]
    stacked <- .VStack(args = lapply(diffs, function(diff) { Reshape(diff, rows = 1, cols = length) }))
    SumEntries(Norm(stacked, p = "fro", axis = 2))
  }
}
