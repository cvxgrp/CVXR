#'
#' The Atom class.
#'
#' This virtual class represents atomic expressions in CVXR.
#'
#' @aliases Atom
#' @export
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

setMethod("validate_args", "Atom", function(object) { })
setMethod("size", "Atom", function(object) { object@.size })
setMethod("is_positive", "Atom", function(object) { sign_from_args(object)[1] })
setMethod("is_negative", "Atom", function(object) { sign_from_args(object)[2] })
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

setMethod("variables", "Atom", function(object) {
  var_list <- lapply(object@args, function(arg) { variables(arg) })
  unique(flatten_list(var_list))
})

setMethod("parameters", "Atom", function(object) {
  param_list <- lapply(object@args, function(arg) { parameters(arg) })
  unique(flatten_list(param_list))
})

setMethod("constants", "Atom", function(object) {
  const_list <- lapply(object@args, function(arg) { constants(arg) })
  unique(flatten_list(const_list))   # TODO: Is this the correct way to remove duplicates?
})

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
      if(is.na(arg_val) && !is_constant(object))
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

.grad.Atom <- function(object) { stop("Unimplemented") }
setMethod("grad", "Atom", function(object) {
  # Short-circuit to all zeros if known to be constant
  if(is_constant(object))
    return(constant_grad(object))
  
  # Returns NA if variable values are not supplied
  arg_values <- list()
  for(arg in object@args) {
    arg_val <- value(arg)
    if(is.na(arg_val))
      return(error_grad(object))
    else
      arg_values <- c(arg_values, arg_val)
  }
  
  # A list of gradients wrt arguments
  grad_self <- .grad(arg_values)
  
  # The chain rule
  result <- list()
  idx <- 1
  for(arg in object@args) {
    # A dictionary of gradients wrt variables
    # Partial argument / partial x
    grad_arg <- grad(arg)
    for(key in names(grad_arg)) {   # TODO: Check if this usage of key matches CVXPY
      # None indicates gradient is not defined
      if(is.na(grad_arg[key]) || is.na(grad_self[idx]))
        result[key] <- NA
      else {
        D <- grad_arg[key] * grad_self[idx]
        # Convert 1x1 matrices to scalars
        if((is.matrix(D) || is(D, "Matrix")) && dim(D) == c(1,1))
          D <- D[1,1]
        
        if(key %in% names(result))
          result[key] <- result[key] + D
        else
          result[key] <- D
      }
    }
  }
  return(result)
})

.domain.Atom <- function(object) { list() }
setMethod("domain", "Atom", function(object) {
  doms <- lapply(object@args, function(arg) { domain(arg) })
  cons <- lapply(doms, function(dom) { dom })
  c(.domain(object), cons)
})

#'
#' The AxisAtom class.
#'
#' This virtual class represents atomic expressions that can be applied along an axis in CVXR.
#'
#' @slot expr A numeric element, data.frame, matrix, vector, or Expression.
#' @slot axis An integer specifying the axis across which to apply the atom. For a matrix, 1 indicates rows, 2 indicates columns, and NA indicates rows and columns (all elements).
#' @aliases AxisAtom
#' @export
AxisAtom <- setClass("AxisAtom", representation(expr = "ConstValORExpr", axis = "numeric"), prototype(axis = NA_real_), contains = c("VIRTUAL", "Atom"))

setMethod("initialize", "AxisAtom", function(.Object, ..., expr, axis) {
  .Object@expr <- expr
  .Object@axis <- axis
  .Object <- callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("size_from_args", "AxisAtom", function(object) {
  if(is.na(object@axis))
    c(1, 1)
  else if(object@axis == 1)
    c(size(object@args[[1]])[1], 1)
  else   # axis == 2
    c(1, size(object@args[[1]])[2])
})

setMethod("get_data", "AxisAtom", function(object) { list(object@axis) })

setMethod("validate_args", "AxisAtom", function(object) {
  if(length(object@axis) != 1 || !(is.na(object@axis) || object@axis %in% c(1, 2)))
     stop("Invalid argument for axis")
})

.axis_grad.AxisAtom <- function(object, values) {
  m <- size(object@args[[1]])[1]
  n <- size(object@args[[1]])[2]
  if(is.na(object@axis)) {
    value <- matrix(t(values[[1]]), nrow = m*n, ncol = 1)
    D <- .column_grad(object, value)
    if(!is.na(D))
      D <- Matrix(D, sparse = TRUE)
  } else {
    if(object@axis == 2) {   # Function apply to each column
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, n))
      for(i in 1:n) {
        value <- values[[1]][,i]
        d <- t(.column_grad(object, value))
        if(is.na(d))
          return(list(NA))
        row <- seq(i*n, i*n+m-1)
        col <- rep(1,m) * i
        D <- D + sparseMatrix(i = row, j = col, x = as.vector(d), dims = c(m*n, n))
      }
    } else {   # Function apply to each row
      values <- t(values[[1]])
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, m))
      for(i in 1:m) {
        value <- values[,i]
        d <- t(.column_grad(value))
        if(is.na(d))
          return(list(NA))
        row <- seq(i, i+(n-1)*m)
        col <- rep(1,n)*i
        D <- D + sparseMatrix(i = row, j = col, x = as.vector(d), dims = c(m*n, m))
      }
    }
  }
  list(D)
}
.column_grad.AxisAtom <- function(object, value) { stop("Unimplemented") }

.AffineProd <- setClass("AffineProd", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")
AffineProd <- function(x, y) { .AffineProd(x = x, y = y) }

setMethod("initialize", "AffineProd", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., args = list(x, y))
})

setMethod("validate_args", "AffineProd", function(object) {
  if(!is_affine(object@args[[1]]) || !is_affine(object@args[[2]]))
    stop("The arguments to AffineProd must be affine")
  mul_shapes(size(object@args[[1]]), size(object@args[[2]]))
})

setMethod("to_numeric", "AffineProd", function(object) { values[[1]] %*% values[[2]] })
setMethod("size_from_args", "AffineProd", function(object) { mul_shapes(size(object@args[[1]]), size(object@args[[2]])) })
setMethod("sign_from_args", "AffineProd", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })
setMethod("is_atom_convex", "AffineProd", function(object) { FALSE })
setMethod("is_atom_concave", "AffineProd", function(object) { FALSE })
setMethod("is_incr", "AffineProd", function(object, idx) { is_positive(object@args[[2-idx]]) })
setMethod("is_decr", "AffineProd", function(object, idx) { is_negative(object@args[[2-idx]]) })
setMethod("is_quadratic", "AffineProd", function(object) { TRUE })

.grad.AffineProd <- function(object, values) {
  X <- values[[1]]
  Y <- values[[2]]
  
  DX_rows <- prod(size(object@args[[1]]))
  cols <- size(object@args[[1]])[1] * size(object@args[[2]])[2]
  
  # TODO: Finish AffineProd gradient implementation
  list(DX, DY)
}

.GeoMean <- setClass("GeoMean", representation(x = "Expression", p = "numeric", max_denom = "numeric"),
                                prototype(p = NA_real_, max_denom = 1024), contains = "Atom")
GeoMean <- function(x, p = NA_real_, max_denom = 1024) { .GeoMean(x = x, p = p, max_denom  = max_denom) }

# TODO: Finish implementing GeoMean. Need to handle fractions properly and add slots for tree, cone_lb, etc
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
  
  if(is.na(p))
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

setMethod("validate_args", "GeoMean", function(object) { })
setMethod("to_numeric", "GeoMean", function(object, values) {
  values <- as.vector(values[[1]])
  val <- 1.0
  for(idx in length(values)) {
    x <- values[[idx]]
    p <- object@w[idx]
    val <- val * x^p
  }
  val
})

.domain.GeoMean <- function(object) { list(object@args[[1]][object@w > 0] >= 0) }
.grad.GeoMean <- function(object, values) {
  x <- as.matrix(values[[1]])
  # No special case when only one non-zero weight
  w_arr <- as.vector(object@w)
  # Outside domain
  if(any(x[w_arr > 0] <= 0))
    return(list(NA))
  else {
    D <- w_arr/as.vector(x) * to_numeric(object, values)
    return(list(t(Matrix(D, sparse = TRUE))))
  }
}

setMethod("size_from_args", "GeoMean", function(object) { c(1,1) })
setMethod("sign_from_args", "GeoMean", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "GeoMean", function(object) { FALSE })
setMethod("is_atom_concave", "GeoMean", function(object) { TRUE })
setMethod("is_incr", "GeoMean", function(object, idx) { TRUE })
setMethod("is_decr", "GeoMean", function(object, idx) { FALSE })
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

GeoMean.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  w <- data[[1]]
  w_dyad <- data[[2]]
  tree <- data[[3]]
  t <- create_var(c(1,1))
  
  if(size(arg_objs[[1]])[2] == 1)
    x_list <- lapply(1:length(w), function(i) { Index.get_index(arg_objs[[1]], list(), i, 1) })
  else if(size(arg_objs[[1]])[1] == 1)
    x_list <- lapply(1:length(w), function(i) { Index.get_index(arg_objs[[1]], list(), 1, i) })
  list(t, gm_constrs(t, x_list, w))
}

setMethod("graph_implementation", "GeoMean", function(object, arg_objs, size, data = NA_real_) {
  GeoMean.graph_implementation(arg_objs, size, data)
})

HarmonicMean <- function(x) {
  x <- as.Constant(x)
  prod(size(x)) * Pnorm(x = x, p = -1)
}

.LambdaMax <- setClass("LambdaMax", representation(A = "ConstValORExpr"), contains = "Atom")
LambdaMax <- function(A) { .LambdaMax(A = A) }

setMethod("initialize", "LambdaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., args = list(.Object@A))
})

setMethod("validate_args", "LambdaMax", function(object) {
  if(size(object@args[[1]])[1] != size(object@args[[1]])[2])
    stop("The argument to LambdaMax must resolve to a square matrix")
})

setMethod("to_numeric", "LambdaMax", function(object, values) {
  if(!all(t(values[[1]]) == values[[1]]))
    stop("LambdaMax called on a non-symmetric matrix")
  max(eigen(values[[1]], only.values = TRUE)$values)
})

setMethod("size_from_args", "LambdaMax", function(object) { c(1, 1) })
setMethod("sign_from_args", "LambdaMax", function(object) { c(FALSE, FALSE) })
setMethod("is_atom_convex", "LambdaMax", function(object) { FALSE })
setMethod("is_atom_concave", "LambdaMax", function(object) { FALSE })
setMethod("is_incr", "LambdaMax", function(object, idx) { FALSE })
setMethod("is_decr", "LambdaMax", function(object, idx) { FALSE })

.domain.LambdaMax <- function(object) { list(t(object@args[[1]]) == object@args[[1]]) }
.grad.LambdaMax <- function(object, values) {
  r <- eigen(values[[1]], only.values = FALSE)
  v <- r$vectors  # eigenvectors
  w <- r$values   # eigenvalues
  
  d <- rep(0, length(w))
  d[length(d)] <- 1
  d <- diag(d)
  D <- v %*% d %*% t(v)
  list(t(Matrix(as.vector(D), sparse = TRUE)))
}

LambdaMax.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]
  n <- size(A)[1]
  # SDP constraint.
  t <- create_var(c(1,1))
  prom_t <- promote(t, c(n,1))
  # I*t - A
  expr <- sub_expr(diag_vec(prom_t), A)
  list(t, list(SDP(expr)))
}

setMethod("graph_implementation", "LambdaMax", function(object, arg_objs, size, data = NA_real_) {
  LambdaMax.graph_implementation(arg_objs, size, data)
})

LambdaMin <- function(X) {
  X <- as.Constant(X)
  -LambdaMax(-X)
}

LambdaSumLargest <- function(X, k) {
  X <- as.Constant(X)
  if(size(X)[1] != size(X)[2])
    stop("First argument must be a square matrix")
  else if(as.integer(k) != k || k <= 0)
    stop("Second argument must be a positive integer")
  
  Z <- Semidef(size(X)[1])
  k*LambdaMax(X - Z) + Trace(Z)
}

LambdaSumSmallest <- function(X, k) {
  X <- as.Constant(X)
  -LambdaSumLargest(-X, k)
}

.LogDet <- setClass("LogDet", representation(A = "ConstValORExpr"), contains = "Atom")
LogDet <- function(A) { .LogDet(A = A) }

setMethod("initialize", "LogDet", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., args = list(.Object@A))
})

setMethod("validate_args", "LogDet", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("The argument to LogDet must be a square matrix")
})

setMethod("to_numeric", "LogDet", function(object, values) {
  logdet <- determinant(values[[1]], logarithm = TRUE)
  if(logdet$sign == 1)
    return(as.numeric(logdet$modulus))
  else
    return(-Inf)
})

setMethod("size_from_args", "LogDet", function(object) { c(1, 1) })
setMethod("sign_from_args",  "LogDet", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "LogDet", function(object) { FALSE })
setMethod("is_atom_concave", "LogDet", function(object) { TRUE })
setMethod("is_incr", "LogDet", function(object, idx) { FALSE })
setMethod("is_decr", "LogDet", function(object, idx) { FALSE })

.grad.LogDet <- function(object, values) {
  X <- as.matrix(values[[1]])
  eigen_val <- eigen(X, only.values = TRUE)$values
  if(min(eigen_val) > 0) {
    # Grad: t(X^(-1))
    D <- t(solve(X))
    return(list(t(Matrix(as.vector(D), sparse = TRUE))))
  } else   # Outside domain
    return(list(NA))
}
.domain.LogDet <- function(object) { list(object@args[[1]] %>>% 0) }

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
  constraints <- c(constraints, SDP(A))
  
  # Fix Z as upper triangular, D as diagonal, and diag(D) as diag(Z)
  Z_lower_tri <- upper_tri(transpose(Z))
  constraints <- c(constraints, create_eq(Z_lower_tri))
  
  # D[i,i] = Z[i,i]
  constraints <- c(constraints, create_eq(D, diag_mat(Z)))
  
  # Fix X using the fact that A must be affine by the DCP rules
  # X[1:n, 1:n] == D
  constraints <- Index.block_eq(X, diag_vec(D), constraints, 1, n, 1, n)
  
  # X[1:n, n:2*n] == Z
  constraints <- Index.block_eq(X, Z, constraints, 1, n, n, 2*n)
  
  # X[n:2*n, n:2*n] == A
  constraints <- Index.block_eq(X, A, constraints, n, 2*n, n, 2*n)
  
  # Add the objective sum(log(D[i,i]))
  graph <- Log.graph_implementation(list(D), c(n, 1))
  obj <- graph[[1]]
  constr <- graph[[2]]
  list(sum_entries(obj), c(constraints, constr))
}

setMethod("graph_implementation", "LogDet", function(object, arg_objs, size, data = NA_real_) {
  LogDet.graph_implementation(arg_objs, size, data)
})

.LogSumExp <- setClass("LogSumExp", contains = "AxisAtom")
LogSumExp <- function(x, axis = NA_real_) { .LogSumExp(expr = x, axis = axis) }

setMethod("to_numeric", "LogSumExp", function(object, values) {
  if(is.na(object@axis))
    log(sum(exp(values[[1]])))
  else
    log(apply(exp(values[[1]]), object@axis, sum))
})

.grad.LogSumExp <- function(object, values) { .axis_grad(object, values) }
.column_grad.LogSumExp <- function(object, value) {
  denom <- log(sum(exp(value)))
  nom <- exp(value)
  D <- nom/denom
  D
}

setMethod("sign_from_args",  "LogSumExp", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "LogSumExp", function(object) { TRUE })
setMethod("is_atom_concave", "LogSumExp", function(object) { FALSE })
setMethod("is_incr", "LogSumExp", function(object, idx) { FALSE })
setMethod("is_decr", "LogSumExp", function(object, idx) { FALSE })

LogSumExp.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  axis <- data[[1]]
  t <- create_var(size)
  
  # sum(exp(x - t)) <= 1
  if(is.na(axis)) {
    prom_t <- promot(t, size(x))
    expr <- sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]
    obj <- sum_entries(obj)
  } else if(axis == 2) {
    prom_size <- c(size(x)[1], 1)
    ones <- create_const(matrix(1, nrow = prom_size[1], ncol = prom_size[2]), prom_size)
    prom_t <- mul_expr(ones, t, size(x))
    expr <- sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]
    
    const_size <- c(1, size(x)[1])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- mul_expr(ones, obj, size)
  } else {    # axis == 1
    prom_size <- c(1, size(x)[2])
    ones <- create_const(matrix(1, nrow = prom_size[1], ncol = prom_size[2]), prom_size)
    prom_t <- rmul_expr(t, ones, size(x))
    expr <- sub_expr(x, prom_t)
    graph <- Exp.graph_implementation(list(expr), size(x))
    obj <- graph[[1]]
    constraints <- graph[[2]]
    
    const_size <- c(size(x)[2], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- rmul_expr(obj, ones, size)
  }

  ones <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  constraints <- c(constraints, create_leq(obj, ones))
  list(t, constraints)
}

setMethod("graph_implementation", "LogSumExp", function(object, arg_objs, size, data = NA_real_) {
  LogSumExp.graph_implementation(arg_objs, size, data)
})

.MatrixFrac <- setClass("MatrixFrac", representation(X = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")
MatrixFrac <- function(X, P) { .MatrixFrac(X = X, P = P) }

setMethod("initialize", "MatrixFrac", function(.Object, ..., X, P) {
  .Object@X <- X
  .Object@P <- P
  callNextMethod(.Object, ..., args = list(.Object@X, .Object@P))
})

setMethod("validate_args", "MatrixFrac", function(object) {
  X <- object@args[[1]]
  P <- object@args[[2]]
  if(size(P)[1] != size(P)[2])
    stop("The second argument to MatrixFrac must be a square matrix")
  else if(size(X)[1] != size(P)[1])
    stop("The arguments to MatrixFrac have incompatible dimensions")
})

setMethod("to_numeric", "MatrixFrac", function(object, values) {
  # TODO: Raise error if not invertible?
  X <- values[[1]]
  P <- values[[2]]
  sum(diag(t(X) %*% solve(P) %*% X))
})

setMethod("size_from_args", "MatrixFrac", function(object) { c(1, 1) })
setMethod("sign_from_args", "MatrixFrac", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "MatrixFrac", function(object) { TRUE })
setMethod("is_atom_concave", "MatrixFrac", function(object) { FALSE })
setMethod("is_incr", "MatrixFrac", function(object, idx) { FALSE })
setMethod("is_decr", "MatrixFrac", function(object, idx) { FALSE })
setMethod("is_quadratic", "MatrixFrac", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

.domain.MatrixFrac <- function(object) { list(object@args[[2]] %>>% 0) }
.grad.MatrixFrac <- function(object, values) {
  X <- as.matrix(values[[1]])
  P <- as.matrix(values[[2]])
  P_inv <- solve(P)
  
  # partial_X = (P^-1+P^-T)X
  # partial_P = (P^-1 * X * X^T * P^-1)^T
  DX <- (P_inv + t(P_inv)) %*% X
  DX <- as.vector(t(DX))
  DX <- t(Matrix(DX, sparse = TRUE))
  
  DP <- P_inv %*% X
  DP <- DP %*% t(X)
  DP <- DP %*% P_inv
  DP <- -t(DP)
  DP <- t(Matrix(as.vector(t(DP))))
  list(DX, DP)
}

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
  
  # M[1:n, n:n+m] == X
  constraints <- Index.block_eq(M, X, constraints, 1, n, n, n+m)
  
  # M[n:n+m, n:n+m] == Tmat
  constraints <- Index.block_eq(M, Tmat, constraints, n, n+m, n, n+m)
  
  # Add SDP constraints.
  list(trace(Tmat), c(constraints, list(SDP(M))))
}

setMethod("graph_implementation", "MatrixFrac", function(object, arg_objs, size, data = NA_real_) {
  MatrixFrac.graph_implementation(arg_objs, size, data)
})

.MaxEntries <- setClass("MaxEntries", contains = "AxisAtom")
MaxEntries <- function(x, axis = NA_real_) { .MaxEntries(expr = x, axis = axis) }

setMethod("to_numeric", "MaxEntries", function(object, values) {
  if(is.na(object@axis))
    max(values[[1]])
  else
    apply(values[[1]], axis, max)
})

setMethod("sign_from_args",  "MaxEntries", function(object) { c(is_positive(object@args[[1]]), is_negative(object@args[[1]])) })
setMethod("is_atom_convex", "MaxEntries", function(object) { TRUE })
setMethod("is_atom_concave", "MaxEntries", function(object) { FALSE })
setMethod("is_incr", "MaxEntries", function(object, idx) { TRUE })
setMethod("is_decr", "MaxEntries", function(object, idx) { FALSE })

.grad.MaxEntries <- function(object, values) { .axis_grad(object, values) }
.column_grad.MaxEntries <- function(object, value) {
  # Grad: 1 for a largest index
  value <- as.vector(value)
  idx <- which.max(value)
  D <- rep(0, length(value))
  D[idx] <- 1
  D
}

MaxEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  axis <- data[[1]]
  if(is.na(axis)) {
    t <- create_var(c(1,1))
    promoted_t <- promote(t, size(arg_objs[[1]]))
  } else if(axis == 2) {
    t <- create_var(c(1, size(arg_objs[[1]])[2]))
    const_size <- c(size(arg_objs[[1]])[1], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    promoted_t <- mul_expr(ones, t, size(arg_objs[[1]]))
  } else {   # axis == 1
    t <- create_var(c(size(arg_objs[[1]])[1], 1))
    const_size <- c(1, size(arg_objs[[1]])[2])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    promoted_t <- rmul_expr(t, ones, size(arg_objs[[1]]))
  }

  constraints <- list(create_leq(arg_objs[[1]], promoted_t))
  list(t, constraints)
}

setMethod("graph_implementation", "MaxEntries", function(object, arg_objs, size, data = NA_real_) {
  MaxEntries.graph_implementation(arg_objs, size, data)
})

MinEntries <- function(x, axis = NA_real_) {
  x <- as.Constant(x)
  -MaxEntries(-x, axis = axis)
}

max.Expression <- function(..., na.rm = FALSE) {
  if(!na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  max_args <- lapply(vals[is_expr], function(expr) { MaxEntries(expr = expr) })
  if(!all(is_expr)) {
    max_num <- max(sapply(vals[!is_expr], function(v) { max(v, na.rm = na.rm) }))
    max_args <- c(max_args, max_num)
  }
  .MaxElemwise(args = max_args)
}

min.Expression <- function(..., na.rm = FALSE) {
  if(!na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  min_args <- lapply(vals[is_expr], function(expr) { MinEntries(expr = expr) })
  if(!all(is_expr)) {
    min_num <- min(sapply(vals[!is_expr], function(v) { min(v, na.rm = na.rm) }))
    min_args <- c(min_args, min_num)
  }
  min_args <- lapply(min_args, function(arg) { -as.Constant(arg) })
  -.MaxElemwise(args = min_args)
}

#'
#' The Pnorm class.
#'
#' This class represents the p-norm expression.
#'
#' @aliases Pnorm
#' @export
.Pnorm <- setClass("Pnorm", representation(p = "numeric", max_denom = "numeric", .approx_error = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_), contains = "AxisAtom")

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

setMethod("validate_args", "Pnorm", function(object) {
  callNextMethod()
  if(!is.na(object@axis) && object@p != 2)
    stop("The axis parameter is only supported for p = 2")
})

setMethod("name", "Pnorm", function(object) { 
  sprintf("%s(%s, %s)", class(object), name(object@args[1]), object@p) 
})

p_norm <- function(x, p) {
  if(p == Inf)
    max(abs(x))
  else if(p == 0)
    sum(x != 0)
  else if(p >= 1)
    sum(abs(x)^p)^(1/p)
  else
    sum(x^p)^(1/p)
}

setMethod("to_numeric", "Pnorm", function(object, values) {
  if(is.na(object@axis))
    values <- as.vector(values[[1]])
  else
    values <- as.matrix(values[[1]])
  
  if(object@p < 1 && any(values < 0))
    return(-Inf)
  
  if(object@p < 0 && any(values == 0))
    return(0)
  
  if(is.na(object@axis))
    retval <- p_norm(values[[1]], object@p)
  else
    retval <- apply(values, object@axis, function(x) { p_norm(x, object@p) })
  retval
})

setMethod("sign_from_args",  "Pnorm", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "Pnorm", function(object) { object@p >= 1})
setMethod("is_atom_concave", "Pnorm", function(object) { object@p < 1 })
setMethod("is_incr", "Pnorm", function(object, idx) { object@p < 1 || (object@p >= 1 && is_positive(object@args[[1]])) })
setMethod("is_decr", "Pnorm", function(object, idx) { object@p >= 1 && is_negative(object@args[[1]]) })
setMethod("get_data", "Pnorm", function(object) { list(object@p, object@axis) })

.grad.Pnorm <- function(object, values) { .axis_grad(values) }
.column_grad.Pnorm <-function(object, value) {
  rows <- prod(size(object@args[[1]]))
  value <- as.matrix(value)
  
  # Outside domain
  if(object@p < 1 && any(value <= 0))
    return(NA)
  D_null <- sparseMatrix(i = c(), j = c(), dims = c(rows, 1))
  if(object@p == 1) {
    D_null <- D_null + (value > 0)
    D_null <- D_null - (value < 0)
    return(t(Matrix(as.vector(D_null), sparse = TRUE)))   # TODO: Is this redundant? Check against CVXPY
  }
  denominator <- p_norm(value, object@p)
  denominator <- denominator^(object@p - 1)
  
  # Subgrad is 0 when denom is 0 (or undefined)
  if(denominator == 0) {
    if(object@p >= 1)
      return(D_null)
    else
      return(NA)
  } else {
    nominator <- value^(object@p - 1)
    frac <- nominator / denominator
    return(matrix(as.vector(frac)))
  }
}

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
      return(list(t, list(SOCAxis(reshape(t, c(prod(t$size), 1)), x, axis))))
    }
  }
  
  if(p == Inf) {
    t_ <- promote(t, x$size)
    return(list(t, list(create_leq(x, t_), create_geq(sum_expr(list(x, t_))))))
  }
  
  # We need absolute value constraint for symmetric convex branches (p >= 1)
  # We alias |x| as x from this point forward to make code pretty
  if(p >= 1) {
    absx <- create_var(x$size)
    constraints <- c(constraints, create_leq(x, absx), create_geq(sum_expr((list(x, absx)))))
    x <- absx
  }
  
  if(p == 1)
    return(list(sum_entries(x), constraints))
  
  # Now take care of remaining convex/concave branches
  # To create rational powers, need new variable r and constraint sum(r) == t
  r <- create_var(x$size)
  t_ <- promote(t, x$size)
  constraints <- c(constraints, create_eq(sum_entries(r), t))
  
  # TODO: Make p a fraction so input weight to gm_constrs is nice tuple of fractions
  if(p < 0)
    constraints <- c(constraints, gm_constrs(t_, list(x, r), c(-p/(1-p), 1/(1-p)) ))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, t_), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, t_)), c(1/p, 1-1/p))
  
  list(t, constraints)
}

setMethod("graph_implementation", "Pnorm", function(object, arg_objs, size, data = NA_real_) {
  Pnorm.graph_implementation(arg_objs, size, data)
})

setMethod("norm", signature(x = "Expression", type = "character"), function(x, type = c("O", "I", "F", "M", "2")) {
  x <- as.Constant(x)
  
  # Norms for scalars same as absolute value
  if(type == "O" || type == "o" || type == "1")   # Maximum absolute column sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 2))
  else if(type == "I" || type == "i")             # Maximum absolute row sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 1))
  else if(type == "F" || type == "f")             # Frobenius norm (Euclidean norm if x is treated as a vector)
    Pnorm(x = x, p = 2, axis = NA_real_)
  else if(type == "M" || type == "m")             # Maximum modulus (absolute value) of all elements in x
    MaxEntries(Abs(x = x))
  else if(type == "2")                            # Spectral norm (largest singular value of x)
    # Sqrt(LambdaMax(A = t(x) %*% x))
    stop("Spectral norm is currently unimplemented")
  else
    stop("Unrecognized type ", type)
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
      SigmaMax(x = x)
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
  vecnorms <- lapply(1:size(X)[1], function(i) { norm(X[i,], p) })
  
  # Outer norms
  Norm(.HStack(args = vecnorms), q)
}

.NormNuc <- setClass("NormNuc", representation(A = "Expression"), contains = "Atom")
NormNuc <- function(A) { .NormNuc(A = A) }

setMethod("initialize", "NormNuc", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., args = list(.Object@A))
})

setMethod("to_numeric", "NormNuc", function(object, values) {
  # Returns the nuclear norm (i.e. the sum of the singular values) of A
  sum(svd(values[[1]])$d)
})

setMethod("size_from_args", "NormNuc", function(object) { c(1, 1) })
setMethod("sign_from_args",  "NormNuc", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "NormNuc", function(object) { TRUE })
setMethod("is_atom_concave", "NormNuc", function(object) { FALSE })
setMethod("is_incr", "NormNuc", function(object, idx) { FALSE })
setMethod("is_decr", "NormNuc", function(object, idx) { FALSE })

.grad.NormNuc <- function(object, values) {
  # Grad: UV^T
  s <- svd(val_svd)
  D <- s$u %*% t(s$v)
  list(t(Matrix(as.vector(D), sparse = TRUE)))
}

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
  # X[1:rows, rows:rows+cols] == A
  constraints <- Index.block_eq(X, A, constraints, 1, rows, rows, rows+cols)
  half <- create_const(0.5, c(1,1))
  trace <- mul_expr(half, trace(X), c(1, 1))
  
  # Add SDP constraint.
  list(trace, c(list(SDP(X)), constraints))
}

setMethod("graph_implementation", "NormNuc", function(object, arg_objs, size, data = NA_real_) {
  NormNuc.graph_implementation(arg_objs, size, data)
})

.decomp_quad <- function(P, cond = NA, rcond = NA) {
  eig <- eigen(P, only.values = FALSE)
  w <- eig$values
  V <- eig$vectors
  
  if(is.na(rcond))
    cond <- rcond
  if(cond %in% c(NA, -1))
    cond <- 1e6 * .Machine$double.eps   # TODO: Check this is doing the correct thing
  
  scale <- max(abs(w))
  w_scaled <- w / scale
  maskp <- w_scaled > cond
  maskn <- w_scaled < -cond
  
  # TODO: Allow indefinite QuadForm
  if(any(maskp) && any(maskn))
    warning("Forming a non-convex expression QuadForm(x, indefinite)")
  M1 <- V[,maskp] %*% sqrt(w_scaled[maskp])
  M2 <- V[,maskn] %*% sqrt(-w_scaled[maskn])
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
    if(length(M1) > 0)
      ret <- ret + scale * SumSquares(Constant(t(M1)) %*% x)
    else if(length(M2) > 0)
      ret <- ret - scale * SumSquares(Cosntant(t(M2)) %*% x)
    return(ret)
  } else
    stop("At least one argument to QuadForm must be constant")
}

#'
#' The QuadOverLin class.
#'
#' This class represents the sum of squared entries in X divided by scalar y.
#' \sum_{i,j} X_{i,j}^2/y
#'
#' @aliases QuadOverLin
#' @export
.QuadOverLin <- setClass("QuadOverLin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")
QuadOverLin <- function(x, y) { .QuadOverLin(x = x, y = y) }

setMethod("initialize", "QuadOverLin", function(.Object, ..., x = .Object@x, y = .Object@y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y))
})

setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@args[[2]]))
    stop("[QuadOverLin: validation] y must be a scalar")
})

setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(values[[1]]^2) / values[[2]] })
setMethod("size_from_args", "QuadOverLin", function(object) { c(1, 1) })
setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })
setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })
setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_positive(object@args[[idx]]) })
setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_negative(object@args[[idx]])) || (idx == 2) })
setMethod("is_quadratic", "QuadOverLin", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

.domain.QuadOverLin <- function(object) { list(object@args[[2]] >= 0) }
.grad.QuadOverLin <- function(object, values) {
  X <- values[[1]]
  y <- values[[2]]
  if(y <= 0)
    return(list(NA, NA))
  else {
    # DX = 2X/y, Dy = -||X||^2_2/y^2
    Dy <- -sum(X^2)/y^2
    Dy <- Matrix(Dy, sparse = TRUE)
    DX <- 2.0*X/y
    DX <- Matrix(as.vector(DX), sparse = TRUE)
    return(list(DX, Dy))
  }
}

QuadOverLin.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  y <- arg_objs[[2]]   # Known to be a scalar.
  v <- create_var(c(1,1))
  two <- create_const(2, c(1,1))
  constraints <- list(SOC(sum_expr(list(y, v)),
                          list(sub_expr(y, v),
                               mul_expr(two, x, x$size))),
                      create_geq(y))
  list(v, constraints)
}

setMethod("graph_implementation", "QuadOverLin", function(object, arg_objs, size, data = NA_real_) {
  QuadOverLin.graph_implementation(arg_objs, size, data)
})

.SigmaMax <- setClass("SigmaMax", representation(A = "Expression"), contains = "Atom")
SigmaMax <- function(A = A) { .SigmaMax(A = A) }

setMethod("initialize", "SigmaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., args = list(.Object@A))
})

setMethod("to_numeric", "SigmaMax", function(object, values) { norm(values[[1]], type = "2") })
setMethod("size_from_args", "SigmaMax", function(object) { c(1, 1) })
setMethod("sign_from_args",  "SigmaMax", function(object) { c(TRUE, FALSE) })
setMethod("is_atom_convex", "SigmaMax", function(object) { TRUE })
setMethod("is_atom_concave", "SigmaMax", function(object) { FALSE })
setMethod("is_incr", "SigmaMax", function(object, idx) { FALSE })
setMethod("is_decr", "SigmaMax", function(object, idx) { FALSE })

.grad.SigmaMax <- function(object, values) {
  # Grad: U diag(e_1) t(V)
  s <- svd(values[[1]])
  ds <- rep(0, length(s$d))
  ds[1] <- 1
  D <- s$u %*% diag(ds) %*% t(s$v)
  t(Matrix(as.vector(D), sparse = TRUE))
}

SigmaMax.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]   # n by m matrix
  size <- dim(A)
  n <- size[1]
  m <- size[2]
  
  # Create a matrix with Schur complement I*t - (1/t)*t(A)*A
  X <- create_var(c(n+m, n+m))
  t <- create_var(c(1,1))
  constraints <- list()
  
  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:n, 1:n] == I_n*t
  prom_t <- promote(t, c(n,1))
  constraints <- Index.block_eq(X, diag_vec(prom_t), constraints, 1, n, 1, n)
  
  # X[1:n, n:n+m] == A
  constraints <- Index.block_eq(X, A, constraints, 1, n, n, n+m)
  
  # X[n:n+m, n:n+m] == I_m*t
  prom_t <- promote(t, c(m,1))
  constraints <- Index.block_eq(X, diag_vec(prom_t), constraints, n, n+m, n, n+m)
  
  # Add SDP constraint.
  list(t, c(constraints, list(SDP(X))))
}

setMethod("graph_implementation", "SigmaMax", function(object, arg_objs, size, data = NA_real_) {
  SigmaMax.graph_implementation(arg_objs, size, data)
})

.SumLargest <- setClass("SumLargest", representation(x = "Expression", k = "numeric"), 
                       validity = function(object) {
                         if(as.integer(object@k) != object@k || object@k <= 0)
                           stop("[SumLargest: validation] k must be a positive integer")
                         return(TRUE)
                         }, contains = "Atom")

SumLargest <- function(x, k) { .SumLargest(x = x, k = k) }

setMethod("initialize", "SumLargest", function(.Object, ..., x, k) {
  .Object@x <- x
  .Object@k <- k
  callNextMethod(.Object, ..., args = list(.Object@x))
})

setMethod("validate_args",   "SumLargest", function(object) {
  if(as.integer(object@k) != object@k || object@k <= 0)
    stop("[SumLargest: validation] k must be a positive integer")
})

setMethod("to_numeric", "SumLargest", function(object, values) {
  # Return the sum of the k largest entries of the matrix
  value <- as.vector(values[[1]])
  k <- min(object@k, length(value))
  val_sort <- sort(value, decreasing = FALSE)
  sum(val_sort[1:k])
})

setMethod("size_from_args", "SumLargest", function(object) { c(1, 1) })
setMethod("sign_from_args", "SumLargest", function(object) { c(is_positive(object@args[[1]]), is_negative(object@args[[1]])) })
setMethod("is_atom_convex", "SumLargest", function(object) { TRUE })
setMethod("is_atom_concave", "SumLargest", function(object) { FALSE })
setMethod("is_incr", "SumLargest", function(object, idx) { TRUE })
setMethod("is_decr", "SumLargest", function(object, idx) { FALSE })
setMethod("get_data", "SumLargest", function(object) { list(object@k) })

.grad.SumLargest <- function(object, values) {
  # Grad: 1 for each of the k largest indices
  value <- as.vector(values[[1]])
  k <- min(object@k, length(value))
  indices <- order(value, decreasing = TRUE)
  D <- rep(0, prod(size(object@args[[1]])))
  D[indices[1:k]] <- 1
  list(Matrix(D, sparse = TRUE))
}

SumLargest.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # min SumEntries(t) + k*q
  # s.t. x <= t + q, t >= 0
  x <- arg_objs[[1]]
  k <- create_const(data[[1]], c(1,1))
  q <- create_var(c(1,1))
  t <- create_var(x$size)
  
  sum_t <- sum_entries(t)
  obj <- sum_expr(list(sum_t, mul_expr(k, q, c(1,1))))
  prom_q <- promote(q, x$size)
  
  constr <- c(create_leq(x, sum_expr(list(t, prom_q))), create_geq(t))
  list(obj, constr)
}

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
    Norm(value[-1] - value[1:(max(rows, cols)-1)], 1)
  else {   # L2 norm for matrices
    args <- lapply(list(...), function(arg) { as.Constant(arg) })
    values <- c(list(value), args)
    diffs <- list()
    for(mat in values) {
      diffs <- c(diffs, list(mat[1:(rows-1), 2:cols] - mat[1:(rows-1), 1:(cols-1)],
                             mat[2:rows, 1:(cols-1)] - mat[1:(rows-1), 1:(cols-1)]))
    }
    length <- size(diffs[[1]])[1] * size(diffs[[2]])[2]
    stacked <- .VStack(args = lapply(diffs, function(diff) { matrix(diff, nrow = 1, ncol = length) }))
    SumEntries(Pnorm(stacked, p = "fro", axis = 1))
  }
}
