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

setMethod("initialize", "Atom", function(.Object, ..., .args = list()) {
  # excl_names = c(".args", slotNames("Expression"))
  # .Object@.args = get_slots(.Object, exclude = excl_names)
  .Object@.args <- lapply(.args, as.Constant)
  .Object@dcp_attr <- init_dcp_attr(.Object)
  validate_args(.Object)
  callNextMethod(.Object, ...)
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
    arg_objs <- lapply(object@.args, function(arg) { canonical_form(arg)[[1]] })
    constraints <- lapply(object@.args, function(arg) { canonical_form(arg)[[2]] })
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
    result <- matrix(rep(0, size[1] * size[2]), nrow = size[1], ncol = size[2])
  } else {
    arg_values <- list()
    for(arg in object@.args) {
      if(is.na(value(arg)) && !is_constant(arg))
        return(NA)
      else
        arg_values <- c(arg_values, value(arg))
    }
    result <- lapply(arg_values, function(x) { as.matrix(x) })
  }
  
  if(length(result) == 1 && all(dim(result[[1]]) == c(1,1)))
    as.numeric(result[[1]])
  else
    result
})

HarmonicMean <- function(x) {
  x <- as.Constant(x)
  prod(size(x)) * Pnorm(x = x, p = -1)
}

KLDiv <- setClass("KLDiv", representation(x = "Expression", y = "Expression"), contains = "Atom")

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
  t <- create_var(c(1,1))
  prom_t <- promote(t, c(n,1))
  expr <- sub_expr(diag_vec(prom_t), A)
  list(t, list(SDP(expr)))
}

.LambdaMin <- setClass("LambdaMin", representation(A = "ConstValORExpr"), contains = "Atom")
LambdaMin <- function(A) { .LambdaMin(A = A) }

setMethod("validate_args", "LambdaMin", function(object) {
  if(size(object@.args[[1]])[1] != size(object@.args[[2]])[2])
    stop("The argument to LambdaMin must resolve to a square matrix")
})

setMethod("initialize", "LambdaMin", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., .args = list(.Object@A))
})

setMethod("shape_from_args", "LambdaMin", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "LambdaMin", function(object) { Sign.UNKNOWN })
setMethod("func_curvature", "LambdaMin", function(object) { Curvature.CONCAVE })
setMethod("monotonicity", "LambdaMin", function(object) { NONMONOTONIC })

LambdaMin.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  A <- arg_objs[[1]]
  n <- size(A)[1]
  t <- create_var(c(1,1))
  prom_t <- promote(t, c(n,1))
  expr <- sub_expr(A, diag_vec(prom_t))
  list(t, list(SDP(expr)))
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

LogDet <- setClass("LogDet", representation(A = "matrix"), contains = "Atom")

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

# LogDet.graph_implementation <- function(arg_objs, size, data = NA_real_) {
#   A <- arg_objs[[1]]  # n by n matrix
#   n <- size(A)[1]
#   X <- create_var(c(2*n, 2*n))
#   canon <- canonical_form(Semidef(2*n))
#   X <- canon[[1]]
#   constraints <- canon[[2]]
#   Z <- create_var(c(n,n))
#   D <- create_var(c(n,1))
#   
#   # Require that X and A are PSD
#   constraints <- c(constraints, SDP(A))
#   
#   # Fix Z as upper triangular, D as diagonal, and diag(D) as diag(Z)
#   Z_lower_tri <- upper_tri(transpose(Z))
#   constraints <- c(constraints, create_eq(Z_lower_tri))
#   
#   # D[i,i] = Z[i,i]
#   constraints <- c(constraints, create_eq(D, diag_mat(Z)))
#   
#   # Fix X using the fact that A must be affine by the DCP rules
#   # X[1:n, 1:n] == D
#   block_eq(X, diag_vec(D), constraints, 1, n, 1, n)
#   
#   # X[1:n, 1:n] == Z
#   block_eq(X, Z, constraints, 1, n, n, 2*n)
#   
#   # X[n:2*n, n:2*n] == A
#   block_eq(X, A, constraints, n, 2*n, n, 2*n)
# }

LogSumExp <- setClass("LogSumExp", representation(x = "Expression"), contains = "Atom")

setMethod("initialize", "LogSumExp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("shape_from_args", "LogSumExp", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "LogSumExp", function(object) { Sign.UNKNOWN })
setMethod("func_curvature",  "LogSumExp", function(object) { Curvature.CONVEX })
setMethod("monotonicity",    "LogSumExp", function(object) { INCREASING })

LogSumExp.graph_implementation <- function(object) {
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  
  # sum(exp(x - t))
  prom_t <- promote(t, size(x))
  expr <- sub_expr(x, prom_t)
  graph0 <- Exp.graph_implementation(list(expr), size(x))
  obj <- graph0[[1]]
  constraints <- graph0[[2]]
  
  graph1 <- SumEntries.graph_implementation(list(obj), c(1,1))
  obj <- graph1[[1]]
  constr <- graph1[[2]]
  
  # obj <= 1
  one <- create_const(1, c(1,1))
  constraints <- c(constraints, constr, create_leq(obj, one))
  list(t, constraints)
}

MatrixFrac <- setClass("MatrixFrac", representation(x = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

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
  # TODO: Finish once index graph implementation done
  list(t, c(constraints, SDP(M)))
}

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

.MinEntries <- setClass("MinEntries", contains = "MaxEntries")
MinEntries <- function(x) { .MinEntries(x = x) }

setMethod("func_curvature", "MinEntries", function(object) { Curvature.CONCAVE })

MinEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(c(1,1))
  promoted_t <- promote(t, size(x))
  constraints <- list(create_leq(promoted_t, x))
  list(t, constraints)
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
setMethod("func_curvature",  "Pnorm", function(object) { Curvature(curvature = ifelse(object@p >= 1, CURV_CONVEX_KEY, CURV_CONCAVE_KEY)) })
setMethod("monotonicity",    "Pnorm", function(object) { ifelse(object@p >= 1, SIGNED, INCREASING) })

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
    t_ <- promote(t, size(x))
    return(list(t, list(create_leq(x, t_), create_geq(sum_expr(list(x, t_))))))
  }
  
  if(p >= 1) {
    absx <- create_var(size(x))
    constraints <- c(constraints, create_leq(x, absx), create_geq(sum_expr((list(x, absx)))))
    x <- absx
  }
  
  if(p == 1)
    return(list(sum_entries(x), constraints))
  
  r <- create_var(size(x))
  t_ <- promote(t, size(x))
  constraints <- c(constraints, create_eq(sum_entries(r), t))
  
  if(p < 0)
    constraints <- c(constraints, gm_constrs(t_, list(x, r), c(-p/(1-p), 1/(1-p)) ))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, t_), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, t_)), c(1/p, 1-1/p))
  
  list(t, constraints)
}

setMethod("norm", signature(x = "Expression", type = "numeric"), function(x, type) { Pnorm(x = x, p = type) })
Norm1   <- function(x) { Pnorm(x = x, p = 1) }
Norm2   <- function(x) { Pnorm(x = x, p = 2) }
NormInf <- function(x) { Pnorm(x = x, p = Inf) }
# TODO: MixedNorm implementation
NormNuc <- setClass("NormNuc", representation(A = "Expression"), contains = "Atom")

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
  
  X <- create_var(c(rows+cols, rows+cols))
  constraints <- list()
  # TODO: Finish implementation once index graph implementation is done
  list()
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
  y <- arg_objs[[2]]
  v <- create_var(c(1,1))
  two <- create_const(2, c(1,1))
  constraints <- list(SOC(sum_expr(c(y, v)),
                          list(sub_expr(y, v),
                               mul_expr(two, x, size(x)))),
                      create_geq(y))
  list(v, constraints)
}

SigmaMax <- setClass("SigmaMax", representation(A = "Expression"), contains = "Atom")

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
  
  X <- create_var(c(n+m, n+m))
  t <- create_var(c(1,1))
  constraints <- list()
  
  prom_t <- promote(t, c(n,1))
  # TODO: Finish implementation once index graph implementation is done
  list()
}

SumLargest <- setClass("SumLargest", representation(x = "Expression", k = "numeric"), 
                       validity = function(object) {
                         if(round(object@k) != object@k || object@k <= 0)
                           stop("[SumLargest: validation] k must be a positive integer")
                         }, contains = "Atom")

setMethod("shape_from_args", "SumLargest", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "SumLargest", function(object) { object@.args[[1]]@dcp_attr@sign })
setMethod("func_curvature", "SumLargest", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "SumLargest", function(object) { INCREASING })

SumLargest.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  k <- create_const(data[[1]], c(1,1))
  q <- create_var(c(1,1))
  t <- create_var(size(x))
  
  graph <- SumEntries.graph_implementation(list(t), c(1,1))
  sum_t <- graph[[1]]
  constr <- graph[[2]]
  obj <- sum_expr(list(sum_t, mul_expr(k, q, c(1,1))))
  prom_q <- promote(q, size(x))
  
  constr <- c(constr, create_leq(x, sum_expr(list(t, prom_q))))
  constr <- c(constr, create_geq(t))
  list(obj, constr)
}

SumSmallest <- function(x, k) {
  x <- as.Constant(x)
  -SumLargest(x = -x, k = k)
}

SumSquares <- function(expr) { QuadOverLin(x = expr, y = 1) }

TotalVariation <- function(value, ...) {
  value <- as.Constant(value)
  rows <- size(values)[1]
  cols <- size(values)[2]
  if(is_scalar(value))
    stop("TotalVariation cannot take a scalar argument")
  else if(is_vector(value))
    norm(value[-1] - value[1:(max(rows, cols) - 1)], 1)
  else {
    args <- lapply(list(...), function(arg) { as.Constant(arg) })
    values <- c(value, args)
    # TODO: Finish when Norm2Elemwise is done
  }
}
