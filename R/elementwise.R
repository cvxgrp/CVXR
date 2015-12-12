Elementwise <- setClass("Elementwise", contains = c("VIRTUAL", "Atom"))

setMethod("validate_args", "Elementwise", function(object) {
  arg_shapes <- lapply(object@.args, function(arg) { arg@dcp_attr@shape })
  Reduce("+", arg_shapes)
})

setMethod("shape_from_args", "Elementwise", function(object) {
  obj_shapes <- lapply(object@.args, function(x) { x@dcp_attr@shape })
  Reduce("+", obj_shapes)
})

.Abs <- setClass("Abs", representation(x = "Expression"), contains = "Elementwise")
setMethod("abs", "Expression", function(x) { .Abs(x = x) })
setMethod("initialize", "Abs", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Abs", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Abs", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Abs", function(object) { SIGNED })
Abs.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size(x))
  constraints <- list(create_geq(sum_expr(list(x, t))), create_leq(x, t))
  list(t, constraints)
}

.Entr <- setClass("Entr", representation(x = "ConstValORExpr"), contains = "Elementwise")
Entr <- function(x) { .Entr(x = x) }

setMethod("initialize", "Entr", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Entr", function(object) { Sign(sign = SIGN_UNKNOWN_KEY) })
setMethod("func_curvature", "Entr", function(object) { Curvature(curvature = CURV_CONCAVE_KEY) })
setMethod("monotonicity", "Entr", function(object) { NONMONOTONIC })
Entr.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(rep(1, size)), size)
  list(t, list(ExpCone(t, x, ones)))
}

.Exp <- setClass("Exp", representation(x = "Expression"), contains = "Elementwise")
setMethod("exp", "Expression", function(x) { .Exp(x = x) })

setMethod("initialize", "Exp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Exp", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Exp", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Exp", function(object) { INCREASING })
Exp.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(rep(1, size)), size)
  list(t, list(ExpCone(x, ones, t)))
}

Huber <- setClass("Huber", representation(x = "Expression", M = "numeric"), 
                           prototype(M = 1), contains = "Elementwise")

setMethod("validate_args", "Huber", function(object) {
  if(!(is_positive(object@M) && is_constant(object@M) && is_scalar(object@M)))
    stop("M must be a non-negative scalar constant")
})

setMethod("initialize", "Huber", function(.Object, ..., x, M = 1) {
  .Object@M <- as.Constant(M)
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Huber", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Huber", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Huber", function(object) { SIGNED })
Huber.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  M <- data[[1]]
  x <- arg_objs[[1]]
  n <- create_var(size)
  s <- create_var(size)
  two <- create_const(2, c(1,1))
  if(is(M, "Parameter"))
    M <- create_param(M, c(1,1))
  else
    M <- create_const(value(M), c(1,1))
  
  # TODO: Finish this implementation
  list()
}

InvPos <- function(x) { Power(x, -1) }

.Log <- setClass("Log", representation(x = "Expression"), contains = "Elementwise")
setMethod("log", "Expression", function(x) { .Log(x = x) })

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Log", function(object) { Sign(sign = SIGN_UNKNOWN_KEY) })
setMethod("func_curvature", "Log", function(object) { Curvature(curvature = CURV_CONCAVE_KEY) })
setMethod("monotonicity", "Log", function(object) { INCREASING })
Log.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  x <- arg_objs[[1]]
  ones <- create_const(matrix(rep(1, size)), size)
  list(t, list(ExpCone(t, ones, x)))
}

.Log1p <- setClass("Log1p", contains = "Log")
setMethod("log1p", "Expression", function(x) { .Log1p(x = x) })
setMethod("sign_from_args", "Log1p", function(object) { object@.args[[1]]@dcp_attr@sign })
Log1p.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  ones <- create_const(matrix(rep(1, size(x))), size(x))
  xp1 <- sum_expr(list(x, ones))
  Log.graph_implementation(list(xp1), size, data)
}

.Logistic <- setClass("Logistic", representation(x = "Expression"), contains = "Elementwise")
Logistic <- function(x) { .Logistic(x = x) }

setMethod("initialize", "Logistic", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Logistic", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Logistic", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Logistic", function(object) { INCREASING })
Logistic.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  x <- arg_objs[[1]]
  t <- create_var(size)
  
  # log(1 + exp(x)) <= t <=> exp(-t) + exp(x - t) <= 1
  graph0 <- Exp.graph_implementation(list(neg_expr(t)), size)
  obj0 <- graph0[[1]]
  constr0 <- graph0[[2]]
  
  graph1 <- Exp.graph_implementation(list(sub_expr(x, t)), size)
  obj1 <- graph1[[1]]
  constr1 <- graph1[[2]]
  
  lhs <- sum_expr(list(obj0, obj1))
  ones <- create_const(matrix(rep(1, size)), size)
  constr <- c(constr0, constr1, create_leq(lhs, ones))
  list(t, constr)
}

.MaxElemwise <- setClass("MaxElemwise", validity = function(object) {
                           if(is.null(object@.args) || length(object@.args) < 2)
                             stop("[MaxElemwise: validation] args must have at least 2 arguments")
                           return(TRUE)
                         }, contains = "Elementwise")
MaxElemwise <- function(arg1, arg2, ...) { .MaxElemwise(.args = list(arg1, arg2, ...)) }

setMethod("sign_from_args", "MaxElemwise", function(object) {
  arg_signs <- lapply(object@.args, function(arg) { arg@dcp_attr@sign })
  contains <- function(sign, args) {
    if(is.null(args) || length(args) == 0) return(FALSE)
    any(sapply(args, function(arg) { arg == sign }))
  }
  
  if(contains(Sign.POSITIVE, arg_signs))
    max_sign <- Sign.POSITIVE
  else if(contains(Sign.ZERO, arg_signs)) {
    if(contains(Sign.UNKNOWN, arg_signs))
      max_sign <- Sign.POSITIVE
    else
      max_sign <- Sign.ZERO
  } else if(contains(Sign.UNKNOWN, arg_signs))
    max_sign <- Sign.UNKNOWN
  else
    max_sign <- Sign.NEGATIVE
  max_sign
})
setMethod("func_curvature", "MaxElemwise", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "MaxElemwise", function(object) { rep(INCREASING, length(object@.args)) })
MaxElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  constraints <- lapply(arg_objs, function(obj) { 
      if(size(obj) != size) 
        obj <- promote(obj, size)
      create_leq(obj, t) 
    })
  list(t, constraints)
}

.MinElemwise <- setClass("MinElemwise", contains = "MaxElemwise")
MinElemwise <- function(arg1, arg2, ...) { .MinElemwise(.args = list(arg1, arg2, ...)) }

setMethod("sign_from_args", "MinElemwise", function(object) {
  arg_signs <- lapply(object@.args, function(arg) { arg@dcp_attr@sign })
  contains <- function(sign, args) {
    if(is.null(args) || length(args) == 0) return(FALSE)
    any(sapply(args, function(arg) { arg == sign }))
  }
  
  if(contains(Sign.NEGATIVE, arg_signs))
    min_sign <- Sign.NEGATIVE
  else if(contains(Sign.ZERO, arg_signs)) {
    if(contains(Sign.UNKNOWN, arg_signs))
      min_sign <- Sign.NEGATIVE
    else
      min_sign <- Sign.ZERO
  } else if(contains(Sign.UNKNOWN, arg_signs))
    min_sign <- Sign.UNKNOWN
  else
    min_sign <- Sign.POSITIVE
  min_sign
})
setMethod("func_curvature", "MinElemwise", function(object) { Curvature.CONCAVE })
MinElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  constraints <- lapply(arg_objs, function(obj) { 
    if(size(obj) != size) 
      obj <- promote(obj, size)
    create_leq(t, obj)
  })
  list(t, constraints)
}

Neg <- function(x) { -MinElemwise(x, 0) }

Norm2Elemwise <- setClass("Norm2Elemwise", contains = "Elementwise")
setMethod("sign_from_args", "Norm2Elemwise", function(object) { Sign.POSITIVE })
setMethod("func_curvature", "Norm2Elemwise", function(object) { Curvature.CONVEX })
setMethod("monotonicity", "Norm2Elemwise", function(object) { rep(SIGNED, length(object@.args)) })
Norm2Elemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  t <- create_var(size)
  arg_objs <- lapply(arg_objs, function(obj) {
    if(size(obj) != size)
      promote(obj, size)
    else
      obj
  })
  list(t, list(SOCElemwise(t, arg_objs)))
}

Pos <- function(x) { MaxElemwise(x, 0) }

.Power <- setClass("Power", representation(x = "Expression", p = "numeric", max_denom = "numeric", .w = "numeric", .approx_error = "numeric"), 
                          prototype(max_denom = 1024, .w = NA_real_, .approx_error = NA_real_), 
                  validity = function(object) {
                    if(!is.na(object@.w))
                      stop("[Validation: power] .w is an internal variable that should not be set by user")
                    else if(!is.na(object@.approx_error))
                      stop("[Validation: power] .approx_error is an internal variable that should not be set by user")
                    return(TRUE)
                    }, contains = "Elementwise")

Power <- function(x, p, max_denom = 1024) { .Power(x = x, p = p, max_denom = max_denom) }
setMethod("^", signature(e1 = "Expression", e2 = "numeric"), function(e1, e2) { Power(x = e1, p = e2) })

setMethod("initialize", "Power", function(.Object, ..., x, p, max_denom = 1024, .w = NA_real_) {
  p_old <- p
  if(p > 1)
    pw <- pow_high(p, max_denom)
  else if(p > 0 && p < 1)
    pw <- pow_mid(p, max_denom)
  else if(p < 0)
    pw <- pow_neg(p, max_denom)
  
  p <- pw[[1]]
  w <- pw[[2]]
  
  if(p == 1) {
    p <- 1
    w <- NA_real_
  } else if(p == 0) {
    p <- 0
    w <- NA_real_
  }
  
  .Object@p <- as.numeric(p)   # TODO: Need to store this as a fraction object, not a rounded numeric
  .Object@.w <- w
  .Object@.approx_error <- abs(.Object@p - p_old)
  
  .Object@x <- x
  .Object@max_denom <- max_denom
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Power", function(object) {
  if(object@p == 1)
    object@.args[[1]]@dcp_attr@sign
  else
    Sign.POSITIVE
})

setMethod("func_curvature", "Power", function(object) {
  if(object@p == 0)
    Curvature.CONSTANT
  else if(object@p == 1)
    Curvature.AFFINE
  else if(object@p < 0 || object@p > 1)
    Curvature.CONVEX
  else if(object@p > 0 && object@p < 1)
    Curvature.CONCAVE
  else
    Curvature.UNKNOWN
})

setMethod("monotonicity", "Power", function(object) {
  if(object@p ==0)
    INCREASING
  else if(object@p == 1)
    INCREASING
  else if(object@p < 0)
    DECREASING
  else if(object@p > 0 && object@p < 1)
    INCREASING
  else if(object@p > 1) {
    if(is_power2(object@p))
      SIGNED
    else
      INCREASING
  }
  else
    stop("Unknown monotonicity for power p = ", object@p)
})

setMethod("get_data", "Power", function(object) { list(object@p, object@.w) })

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
      
      if(p > 0 && p < 1)
        return(list(t, gm_constrs(t, list(x, one), w)))
      else if(p > 1)
        return(list(t, gm_constrs(x, list(t, one)), w))
      else if(p < 0)
        return(list(t, gm_constrs(one, list(x, t)), w))
      else
        stop("This power is not yet supported")
    }
  }
}

setMethod("graph_implementation", "Power", function(object, arg_objs, size, data = NA_real_) {
  Power.graph_implementation(arg_objs, size, data)
})

Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }

Sqrt <- function(x) { Power(x, 1/2) }
setMethod("sqrt", "Expression", function(x) { Sqrt(x) })

Square <- function(x) { Power(x, 2) }
