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

setMethod("initialize", "Atom", function(.Object, ..., .args = list()) {
  # excl_names = c(".args", slotNames("Expression"))
  # .Object@.args = get_slots(.Object, exclude = excl_names)
  .Object@.args <- lapply(.args, as.Constant)
  .Object@dcp_attr <- init_dcp_attr(.Object)
  callNextMethod(.Object, ...)
})

setMethod("Atom.dcp_curvature", signature(curvature = "Curvature", args = "list", monotonicities = "list"),
          function(curvature, args, monotonicities) {
            if(length(args) != length(monotonicities))
              stop("The number of args must be equal to the number of monotonicities")
            arg_curvatures = mapply(function(arg, monotonicity) {
              dcp_curvature(monotonicity = monotonicity, func_curvature = curvature, arg_sign = arg@dcp_attr@sign, arg_curvature = arg@dcp_attr@curvature)
            }, args, monotonicities)
            Reduce("+", arg_curvatures)
          })

#'
#' The pnorm class.
#'
#' This class represents the p-norm expression in CVXR.
#'
#' @aliases pnorm
#' @export
pnorm <- setClass("pnorm", representation(x = "Variable", p = "numeric", max_denom = "numeric", .approx_error = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_), contains = "Atom")
setMethod("initialize", "pnorm", definition = function(.Object, ..., x, p = 2, max_denom = 1024, .approx_error = NA_real_) {
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

setMethod("shape_from_args", "pnorm", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args", "pnorm", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "pnorm", function(object) { Curvature(curvature = ifelse(object@p >= 1, CURV_CONVEX_KEY, CURV_CONCAVE_KEY)) })
setMethod("monotonicity", "pnorm", function(object) { list(ifelse(object@p >= 1, SIGNED, INCREASING)) })

setMethod("get_data", "pnorm", function(object) { list(p = object@p) })
setMethod("name", "pnorm", function(object) { 
  sprintf("%s(%s, %s)", class(object), name(object@.args[1]), object@p) 
})

pnorm.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  p <- data[1]
  x <- arg_objs[1]
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

norm1 <- function(x) { pnorm(x = x, p = 1) }
norm2 <- function(x) { pnorm(x = x, p = 2) }
normInf <- function(x) { pnorm(x = x, p = Inf) }
setMethod("norm", signature(x = "Expression", type = "numeric"), function(x, type) { pnorm(x = x, p = type) })

#'
#' The quad_over_lin class.
#'
#' This class represents the sum of squared entries in X divided by scalar y in CVXR.
#' \sum_{i,j} X_{i,j}^2/y
#'
#' @aliases quad_over_lin
#' @export
quad_over_lin <- setClass("quad_over_lin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")
setMethod("initialize", "quad_over_lin", function(.Object, ..., x = .Object@x, y = .Object@y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., .args = list(.Object@x, .Object@y))
})

setMethod("shape_from_args", "quad_over_lin", function(object) { Shape(rows = 1, cols = 1) })
setMethod("sign_from_args",  "quad_over_lin", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature",  "quad_over_lin", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity",    "quad_over_lin", function(object) { list(SIGNED, DECREASING) })
setMethod("validate_args",   "quad_over_lin", function(object) {
  if(!is_scalar(object@y))
    stop("[quad_over_lin: validation] y must be a scalar")
})

quad_over_lin.graph_implementation <- function(arg_objs, size, data = NA_real_) {
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

setMethod("sum_squares", "Expression", function(expr) { quad_over_lin(x = expr, y = 1) })
