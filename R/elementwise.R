Elementwise <- setClass("Elementwise", contains = c("VIRTUAL", "Atom"))

setMethod("validate_args", "Elementwise", function(object) {
  tot_shape <- object@.args[[1]]@dcp_attr@shape
  if(length(object@.args) > 1) {
    for(arg in object@.args[[-1]])
      tot_shape <- tot_shape + arg@dcp_attr@shape
  }
})

setMethod("shape_from_args", "Elementwise", function(object) {
  obj_shapes <- lapply(object@.args, function(x) { x@dcp_attr@shape })
  Reduce("+", obj_shapes)
})

Abs <- setClass("Abs", representation(x = "Expression"), contains = "Elementwise")
setMethod("initialize", "Abs", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Abs", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Abs", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Abs", function(object) { SIGNED })

Entr <- setClass("Entr", representation(x = "Expression"), contains = "Elementwise")

setMethod("initialize", "Entr", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Entr", function(object) { Sign(sign = SIGN_UNKNOWN_KEY) })
setMethod("func_curvature", "Entr", function(object) { Curvature(curvature = CURV_CONCAVE_KEY) })
setMethod("monotonicity", "Entr", function(object) { NONMONOTONIC })

Exp <- setClass("Exp", representation(x = "Expression"), contains = "Elementwise")

setMethod("sign_from_args", "Exp", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Exp", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Exp", function(object) { INCREASING })

Huber <- setClass("Huber", representation(x = "Expression", M = "numeric"), 
                           prototype(M = 1), contains = "Elementwise")

setMethod("validate_args", "Huber", function(object) {
  if(!(is_positive(object@M) && is_constant(object@M) && is_scalar(object@M)))
    stop("M must be a non-negative scalar constant")
})

setMethod("initialize", "Huber", function(.Object, ..., x, M = 1) {
  .Object@M <- cast_to_const(M)
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Huber", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Huber", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Huber", function(object) { SIGNED })

InvPos <- function(x) { Power(x, -1) }

Log <- setClass("Log", representation(x = "Expression"), contains = "Elementwise")

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Log", function(object) { Sign(sign = SIGN_UNKNOWN_KEY) })
setMethod("func_curvature", "Log", function(object) { Curvature(curvature = CURV_CONCAVE_KEY) })
setMethod("monotonicity", "Log", function(object) { INCREASING })

Log1p <- setClass("Log1p", contains = "Log")
setMethod("sign_from_args", "Log1p", function(object) { object@.args[[1]]@dcp_attr@sign })

Logistic <- setClass("Logistic", representation(x = "Expression"), contains = "Elementwise")

setMethod("initialize", "Logistic", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

Neg <- function(x) { -MinElemwise(x, 0) }

setMethod("sign_from_args", "Logistic", function(object) { Sign(sign = SIGN_POSITIVE_KEY) })
setMethod("func_curvature", "Logistic", function(object) { Curvature(curvature = CURV_CONVEX_KEY) })
setMethod("monotonicity", "Logistic", function(object) { INCREASING })

Pos <- function(x) { MaxElemwise(x, 0) }

.Power <- setClass("Power", representation(x = "Expression", p = "numeric", max_denom = "numeric", .w = "numeric"), 
                          prototype(max_denom = 1024, .w = NA_real_), 
                  validity = function(object) {
                    if(!is.na(object@.w))
                      stop("[Validation: power] .w is an internal variable that should not be set by user")
                    }, contains = "Elementwise")

Power <- function(x, p, max_denom = 1024) { .Power(x = x, p = p, max_denom = max_denom) }

setMethod("initialize", "Power", function(.Object, ..., x, p, max_denom = 1024, .w = NA_real_) {
  # TODO: Fill in p and w accordingly
  .Object@p <- p
  .Object@x <- x
  .Object@max_denom <- max_denom
  callNextMethod(.Object, ..., .args = list(.Object@x))
})

setMethod("sign_from_args", "Power", function(object) {
  if(object@p == 1)
    object@.args[1]@dcp_attr@sign
  else
    Sign(sign = SIGN_POSITIVE_KEY)
})

setMethod("func_curvature", "Power", function(object) {
  if(object@p == 0)
    Curvature(curvature = CURV_CONSTANT_KEY)
  else if(object@p == 1)
    Curvature(curvature = CURV_AFFINE_KEY)
  else if(object@p < 0 || object@p > 1)
    Curvature(curvature = CURV_CONVEX_KEY)
  else if(object@p > 0 && object@p < 1)
    Curvature(curvature = CURV_CONCAVE_KEY)
  else
    Curvature(curvature = CURV_UNKNOWN_KEY)
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

Scalene <- function(x, alpha, beta) { alpha*Pos(x) + beta*Neg(x) }

Sqrt <- function(x) { Power(x, 1/2) }
setMethod("sqrt", "ConstValORExpr", function(x) { Sqrt(x) })

Square <- function(x) { Power(x, 2) }
