#'
#' The Atom class.
#'
#' This virtual class represents atomic expressions in CVXR.
#'
#' @aliases Atom
#' @export
Atom <- setClass("Atom", representation(.args = "list"), prototype(.args = list()), contains = c("VIRTUAL", "Expression"))

# setGeneric("get_slots", function(object, exclude) { standardGeneric("get_slots") })
# setMethod("get_slots", "Atom", function(object, exclude = c()) {
#   slot_names <- slotNames(object)
#   slot_names <- slot_names[which(!(slot_names %in% exclude))]
#   obj_slots <- lapply(slot_names, function(n) { slot(object, n) })
#   names(obj_slots) <- slot_names
#   obj_slots
# })

setMethod("init_dcp_attr", "Atom", function(object) {
  shape <- shape_from_args(object)
  sign <- sign_from_args(object)
  curvature <- Atom.dcp_curvature(func_curvature(object), object@.args, monotonicity(object))
  DCPAttr(sign = sign, curvature = curvature, shape = shape)
})

setMethod("validate_args", "Atom", function(object) { })

setMethod("initialize", "Atom", function(.Object, ..., dcp_attr, .args = list()) {
  # excl_names = c(".args", slotNames("Expression"))
  # .Object@.args = get_slots(.Object, exclude = excl_names)
  .Object@.args <- lapply(.args, as.Constant)
  .Object@dcp_attr <- init_dcp_attr(.Object)
  validate_args(.Object)
  callNextMethod(.Object, ..., dcp_attr = .Object@dcp_attr)
})

setMethod("Atom.dcp_curvature", signature(curvature = "Curvature", args = "list", monotonicities = "character"),
          function(curvature, args, monotonicities) {
            if(length(args) != length(monotonicities))
              stop("The number of args must be equal to the number of monotonicities")
            arg_curvatures = mapply(function(arg, monotonicity) {
              dcp_curvature(monotonicity = monotonicity, func_curvature = curvature, arg_sign = arg@dcp_attr@sign, arg_curvature = arg@dcp_attr@curvature)
            }, args, monotonicities)
            Reduce("+", arg_curvatures)
          })

setMethod("canonicalize", "Atom", function(object) {
  if(is_constant(object)) {
    if(!is.na(parameters(object)) && length(parameters(object)) > 0) {
      size <- size(object)
      param <- CallbackParam(value(object), size[1], size[2])
      return(canonical_form(param))
    } else
      return(canonical_form(Constant(value(object))))
  } else {
    arg_objs <- list()
    constraints <- list()
    for(arg in object@.args) {
      canon <- canonical_form(arg)
      arg_objs[[length(arg_objs) + 1]] <- canon[[1]]
      constraints <- c(constraints, canon[[2]])
    }
    data <- get_data(object)
    graph <- graph_implementation(object, arg_objs, size(object), data)
    return(list(graph[[1]], c(constraints, graph[[2]])))
  }
})

setMethod("graph_implementation", "Atom", function(object, arg_objs, size, data = NA_real_) {
  stop("Unimplemented")
})

setMethod("variables", "Atom", function(object) {
  var_list <- lapply(object@.args, function(arg) { variables(arg) })
  unique(flatten_list(var_list))
})

setMethod("parameters", "Atom", function(object) {
  param_list <- lapply(object@.args, function(arg) { parameters(arg) })
  unique(flatten_list(param_list))
})

setMethod("value", "Atom", function(object) {
  if(is_zero(object)) {
    size <- size(object)
    result <- matrix(0, nrow = size[1], ncol = size[2])
  } else {
    arg_values <- list()
    idx <- 1
    for(arg in object@.args) {
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
  if(is.list(result) && length(result) == 1 && (length(result[[1]]) == 1 || all(dim(result[[1]]) == c(1,1))))
    result[[1]]
  else
    result
})

AxisAtom <- setClass("AxisAtom", representation(expr = "Expression", axis = "numeric"), prototype(axis = NA_real_), contains = c("VIRTUAL", "Atom"))

setMethod("initialize", "AxisAtom", function(.Object, ..., expr, axis) {
  .Object@axis <- axis
  .Object <- callNextMethod(.Object, ..., .args = list(expr))
})

setMethod("size", "AxisAtom", function(object) {
  if(is.na(object@axis))
    c(1, 1)
  else if(object@axis == 0)
    c(1, size(object@.args[[1]])[2])
  else   # axis == 1
    c(size(object@.args[[1]])[1], 1)
})
setMethod("get_data", "AxisAtom", function(object) { list(object@axis) })

setMethod("validate_args", "AxisAtom", function(object) {
  if(length(object@axis) != 1 || !(is.na(object@axis) || object@axis %in% c(0, 1)))
     stop("Invalid argument for axis")
})

GeoMean <- setClass("GeoMean", representation(x = "Expression", p = "numeric", max_denom = "numeric"),
                               prototype(p = NA_real_, max_denom = 1024), contains = "Atom")

setMethod("initialize", "GeoMean", function(.Object, ..., x, p, max_denom) {
  .Object@x <- x
  .Object <- callNextMethod(.Object, ..., .args = list(.Object@x))
  
  x <- .Object@.args[[1]]
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

setMethod("shape_from_args", "GeoMean", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "GeoMean", function(object) { Sign.POSITIVE })
setMethod("func_curvature", "GeoMean", function(object) { Curvature.CONCAVE })
setMethod("monotonicity", "GeoMean", function(object) { INCREASING })
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

GeoMean.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  w <- data[[1]]
  w_dyad <- data[[2]]
  tree <- data[[3]]
  x_list <- lapply(1:length(w), function(i) { Index.get_index(arg_objs[[1]], list(), i, 1) })
  list(t, gm_constrs(t, x_list, w))
}

setMethod("graph_implementation", "GeoMean", function(object, arg_objs, size, data = NA_real_) {
  GeoMean.graph_implementation(arg_objs, size, data)
})

HarmonicMean <- function(x) {
  x <- as.Constant(x)
  prod(size(x)) * Pnorm(x = x, p = -1)
}

.KLDiv <- setClass("KLDiv", representation(x = "Expression", y = "Expression"), contains = "Atom")
KLDiv <- function(x, y) { .KLDiv(x = x, y = y) }

setMethod("validate_args", "KLDiv", function(object) {
  if(!is_scalar(object@.args[[1]]) || !is_scalar(object@.args[[2]]))
    stop("The arguments to KLDiv must resolve to scalars")
})

setMethod("initialize", "KLDiv", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., .args = list(.Object@x, .Object@y))
})

setMethod("shape_from_args", "KLDiv", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "KLDiv", function(object) { Sign.POSITIVE })
setMethod("func_curvature", "KLDiv", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "KLDiv", function(object) { rep(NONMONOTONIC, length(object@.args)) })

KLDiv.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  y <- arg_objs[[2]]
  t <- create_var(c(1,1))
  constraints <- list(ExpCone(t, x, y), create_geq(y))   # y >= 0
  # -t - x + y
  obj <- sub_expr(y, sum_expr(list(x, y)))
  list(obj, constraints)
}

setMethod("graph_implementation", "KLDiv", function(object, arg_objs, size, data = NA_real_) {
  KLDiv.graph_implementation(arg_objs, size, data)
})

.LambdaMax <- setClass("LambdaMax", representation(A = "ConstValORExpr"), contains = "Atom")
LambdaMax <- function(A) { .LambdaMax(A = A) }

setMethod("validate_args", "LambdaMax", function(object) {
  if(size(object@.args[[1]])[1] != size(object@.args[[1]])[2])
    stop("The argument to LambdaMax must resolve to a square matrix")
})

setMethod("initialize", "LambdaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., .args = list(.Object@A))
})

setMethod("shape_from_args", "LambdaMax", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "LambdaMax", function(object) { Sign.UNKNOWN })
setMethod("func_curvature", "LambdaMax", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "LambdaMax", function(object) { NONMONOTONIC })

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
  else if(round(k) != k || k <= 0)
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

setMethod("validate_args", "LogDet", function(object) {
  n <- size(object@.args[[1]])[1]
  m <- size(object@.args[[1]])[2]
  if(n != m)
    stop("The argument to LogDet must be a square matrix")
})

setMethod("initialize", "LogDet", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., .args = list(.Object@A))
})

setMethod("shape_from_args", "LogDet", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "LogDet", function(object) { Sign.UNKNOWN })
setMethod("func_curvature",  "LogDet", function(object) { Curvature.CONCAVE })
setMethod("monotonicity",    "LogDet", function(object) { NONMONOTONIC })

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
  
  # X[1:n, 1:n] == Z
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

.LogSumExp <- setClass("LogSumExp", representation(x = "Expression"), contains = "Atom")
LogSumExp <- function(x) { .LogSumExp(x = x) }

setMethod("initialize", "LogSumExp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("shape_from_args", "LogSumExp", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "LogSumExp", function(object) { Sign.UNKNOWN })
setMethod("func_curvature",  "LogSumExp", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "LogSumExp", function(object) { INCREASING })

LogSumExp.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  
  # sum(exp(x - t))
  prom_t <- promote(t, size(x))
  expr <- sub_expr(x, prom_t)
  graph <- Exp.graph_implementation(list(expr), size(x))
  obj <- sum_entries(graph[[1]])
  constraints <- graph[[2]]

  # obj <= 1
  one <- create_const(1, c(1,1))
  constraints <- c(constraints, create_leq(obj, one))
  list(t, constraints)
}

setMethod("graph_implementation", "LogSumExp", function(object, arg_objs, size, data = NA_real_) {
  LogSumExp.graph_implementation(arg_objs, size, data)
})

.MatrixFrac <- setClass("MatrixFrac", representation(x = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")
MatrixFrac <- function(x, P) { .MatrixFrac(x = x, P = P) }

setMethod("validate_args", "MatrixFrac", function(object) {
  x <- object@.args[[1]]
  P <- object@.args[[2]]
  if(size(P)[1] != size(P)[2])
    stop("The second argument to MatrixFrac must be a square matrix")
  else if(size(x)[2] != 1)
    stop("The first argument to MatrixFrac must be a column vector")
  else if(size(x)[1] != size(P)[1])
    stop("The arguments to MatrixFrac have incompatible dimensions")
})

setMethod("initialize", "MatrixFrac", function(.Object, ..., x, P) {
  .Object@x <- x
  .Object@P <- P
  callNextMethod(.Object, ..., .args = list(.Object@x, .Object@P))
})

setMethod("shape_from_args", "MatrixFrac", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "MatrixFrac", function(object) { Sign.POSITIVE })
setMethod("func_curvature", "MatrixFrac", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "MatrixFrac", function(object) { rep(NONMONOTONIC, length(object@.args)) })

MatrixFrac.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  P <- arg_objs[[2]]   # n by n matrix
  n <- size(P)[1]
  
  # Create a matrix with Schur complement t - t(x)P^(-1)x
  M <- create_var(c(n+1, n+1))
  t <- create_var(c(1,1))
  constraints <- list()
  
  # Fix M using the fact that P must be affine by the DCP rules.
  # M[1:n, 1:n] == P
  constraints <- Index.block_eq(M, P, constraints, 1, n, 1, n)
  
  # M[1:n, n:n+1] == x
  constraints <- Index.block_eq(M, x, constraints, 1, n, n, n+1)
  
  # M[n:n+1, n:n+1] == t
  constraints <- Index.block_eq(M, t, constraints, n, n+1, n, n+1)
  
  # Add SDP constraints.
  list(t, c(constraints, SDP(M)))
}

setMethod("graph_implementation", "MatrixFrac", function(object, arg_objs, size, data = NA_real_) {
  MatrixFrac.graph_implementation(arg_objs, size, data)
})

.MaxEntries <- setClass("MaxEntries", representation(x = "ConstValORExpr"), contains = "Atom")
MaxEntries <- function(x) { .MaxEntries(x = x) }

setMethod("initialize", "MaxEntries", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("shape_from_args", "MaxEntries", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "MaxEntries", function(object) { object@.args[[1]]@dcp_attr@sign })
setMethod("func_curvature",  "MaxEntries", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "MaxEntries", function(object) { INCREASING })

MaxEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  promoted_t <- promote(t, size(x))
  constraints <- list(create_leq(x, promoted_t))
  list(t, constraints)
}

setMethod("graph_implementation", "MaxEntries", function(object, arg_objs, size, data = NA_real_) {
  MaxEntries.graph_implementation(arg_objs, size, data)
})

MinEntries <- function(x) {
  x <- as.Constant(x)
  -MaxEntries(-x)
}

#'
#' The Pnorm class.
#'
#' This class represents the p-norm expression.
#'
#' @aliases Pnorm
#' @export
.Pnorm <- setClass("Pnorm", representation(x = "ConstValORExpr", p = "numeric", max_denom = "numeric", .approx_error = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_), contains = "Atom")

Pnorm <- function(x, p = 2, max_denom = 1024) { .Pnorm(x = x, p = p, max_denom = max_denom) }

setMethod("initialize", "Pnorm", function(.Object, ..., x, p = 2, max_denom = 1024, .approx_error = NA_real_) {
  .Object@x <- x
  .Object@max_denom <- max_denom
  
  p_old <- p
  .Object@p <- p
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
  #  stop("[pnorm: validation] Invalid value for p ", p)
  
  if(.Object@p == Inf)
    .Object@.approx_error <- 0
  else
    .Object@.approx_error <- abs(.Object@p - p_old)
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("shape_from_args", "Pnorm", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "Pnorm", function(object) { Sign.POSITIVE })
setMethod("func_curvature",  "Pnorm", function(object) { if(object@p >= 1) Curvature.CONVEX else Curvature.CONCAVE })
setMethod("monotonicity",    "Pnorm", function(object) { if(object@p >= 1) SIGNED else INCREASING })
setMethod("get_data", "Pnorm", function(object) { object@p })
setMethod("name", "Pnorm", function(object) { 
  sprintf("%s(%s, %s)", class(object), name(object@.args[1]), object@p) 
})

Pnorm.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  p <- data[[1]]
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  constraints <- list()
  
  if(p == 2)
    return(list(t, list(SOC(t, list(x)))))
  
  if(p == Inf) {
    t_ <- promote(t, x$size)
    return(list(t, list(create_leq(x, t_), create_geq(sum_expr(list(x, t_))))))
  }
  
  if(p >= 1) {
    absx <- create_var(x$size)
    constraints <- c(constraints, create_leq(x, absx), create_geq(sum_expr((list(x, absx)))))
    x <- absx
  }
  
  if(p == 1)
    return(list(sum_entries(x), constraints))
  
  r <- create_var(x$size)
  t_ <- promote(t, x$size)
  constraints <- c(constraints, create_eq(sum_entries(r), t))
  
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

setMethod("norm", signature(x = "Expression", type = "character"), function(x, type = c("O", "I", "F", "N", "2")) {
  x <- as.Constant(x)
  
  if(type == "O" || type == "o" || type == "1")
    return(Pnorm(x = x, p = 1))
  else if(type == "I" || type == "i" || type == "inf")
    return(Pnorm(x = x, p = Inf))
  else if(type == "F" || type == "f" || type == "fro")
    return(Pnorm(x = x, p = 2))
  else if(type == "N" || type == "n" || type == "nuc")
    return(NormNuc(A = x))
  else if(type == "2") {
    if(is.matrix(x))
      return(SigmaMax(x = x))
    else
      return(Pnorm(x = x, p = 2))
  } else
    return(Pnorm(x = x, p = type))
})

Norm1   <- function(x) { Pnorm(x = x, p = 1) }
Norm2   <- function(x) { Pnorm(x = x, p = 2) }
NormInf <- function(x) { Pnorm(x = x, p = Inf) }
MixedNorm <- function(X, p = 2, q = 1) {
  X <- as.Constant(X)
  vecnorms <- lapply(1:size(X)[1], function(i) { norm(X[i,], p) })
  norm(HStack(unlist(vecnorms)), q)
}

.NormNuc <- setClass("NormNuc", representation(A = "Expression"), contains = "Atom")
NormNuc <- function(A) { .NormNuc(A = A) }

setMethod("initialize", "NormNuc", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., .args = list(.Object@A))
})

setMethod("shape_from_args", "NormNuc", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "NormNuc", function(object) { Sign.POSITIVE })
setMethod("func_curvature",  "NormNuc", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "NormNuc", function(object) { NONMONOTONIC })

NormNuc.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]
  rows <- size(A)[1]
  cols <- size(A)[2]
  
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
  list(trace, c(SDP(X), constraints))
}

setMethod("graph_implementation", "NormNuc", function(object, arg_objs, size, data = NA_real_) {
  NormNuc.graph_implementation(arg_objs, size, data)
})

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

setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@.args[[2]]))
    stop("[QuadOverLin: validation] y must be a scalar")
})

setMethod("initialize", "QuadOverLin", function(.Object, ..., x = .Object@x, y = .Object@y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., .args = list(.Object@x, .Object@y))
})

setMethod("shape_from_args", "QuadOverLin", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "QuadOverLin", function(object) { Sign.POSITIVE })
setMethod("func_curvature",  "QuadOverLin", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "QuadOverLin", function(object) { c(SIGNED, DECREASING) })

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
  callNextMethod(.Object, ..., .args = list(.Object@A))
})

setMethod("shape_from_args", "SigmaMax", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "SigmaMax", function(object) { Sign.POSITIVE })
setMethod("func_curvature",  "SigmaMax", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "SigmaMax", function(object) { NONMONOTONIC })

SigmaMax.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]   # n by m matrix
  n <- size(A)[1]
  m <- size(A)[2]
  
  # Create a matrix with Schur complement I*t - (1/t)*t(A)*A
  X <- create_var(c(n+m, n+m))
  t <- create_var(c(1,1))
  constraints <- list()
  
  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:n, 1:n] == I_n*t
  prom_t <- promote(t, c(n,1))
  constraints <- Index.block_eq(X, diag_vec(prom_t), constraints, 1, n, 1, n)
  
  # X[1:n, n:n+1] == A
  constraints <- Index.block_eq(X, A, constraints, 1, n, n, n+m)
  
  # X[n:n+m, n:n+m] == I_m*t
  prom_t <- promote(t, c(m,1))
  constraints <- Index.block_eq(X, diag_vec(prom_t), constraints, n, n+m, n, n+m)
  
  # Add SDP constraint.
  list(t, c(constraints, SDP(X)))
}

setMethod("graph_implementation", "SigmaMax", function(object, arg_objs, size, data = NA_real_) {
  SigmaMax.graph_implementation(arg_objs, size, data)
})

.SumLargest <- setClass("SumLargest", representation(x = "Expression", k = "numeric"), 
                       validity = function(object) {
                         if(round(object@k) != object@k || object@k <= 0)
                           stop("[SumLargest: validation] k must be a positive integer")
                         return(TRUE)
                         }, contains = "Atom")

SumLargest <- function(x, k) { .SumLargest(x = x, k = k) }

setMethod("validate_args",   "SumLargest", function(object) {
  if(round(object@k) != object@k || object@k <= 0)
    stop("[SumLargest: validation] k must be a positive integer")
})

setMethod("initialize", "SumLargest", function(.Object, ..., x, k) {
  .Object@x <- x
  .Object@k <- k
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("shape_from_args", "SumLargest", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "SumLargest", function(object) { object@.args[[1]]@dcp_attr@sign })
setMethod("func_curvature", "SumLargest", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "SumLargest", function(object) { INCREASING })

SumLargest.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  k <- create_const(data[[1]], c(1,1))
  q <- create_var(c(1,1))
  t <- create_var(size(x))
  
  sum_t <- sum_entries(t)
  obj <- sum_expr(list(sum_t, mul_expr(k, q, c(1,1))))
  prom_q <- promote(q, size(x))
  
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
  rows <- size(value)[1]
  cols <- size(value)[2]
  if(is_scalar(value))
    stop("TotalVariation cannot take a scalar argument")
  else if(is_vector(value))   # L1 norm for vectors
    norm(value[-1] - value[1:(max(rows, cols) - 1)], 1)
  else {   # L2 norm for matrices
    args <- lapply(list(...), function(arg) { as.Constant(arg) })
    values <- c(value, args)
    diffs <- lapply(values, function(mat) {
      list(mat[1:(rows-1), 2:cols] - mat[1:(rows-1), 1:(cols-1)],
           mat[2:rows, 1:(cols-1)] - mat[1:(rows-1), 1:(cols-1)])
    })
    SumEntries(Norm2Elemwise(.args = flatten_list(diffs)))
  }
}
