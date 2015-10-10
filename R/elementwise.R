Elementwise <- setClass("Elementwise", contains = "Atom")
setMethod("shape_from_args", "Elementwise", function(object) {
  obj_shapes = lapply(object@args, function(x) { x@dcp_attr@shape })
  Reduce("+", obj_shapes)
})

setMethod("validate_args", "Elementwise", function(object) {
  tot_shape = object@args[[1]]@dcp_attr@shape
  for(arg in object@args[-1])
    tot_shape = tot_shape + arg@dcp_attr@shape
})

power <- setClass("power", representation(x = "Expression", p = "numeric", max_denom = "numeric", .w = "numeric"), 
                          prototype(max_denom = 1024, .w = NA_real_), 
                  validity = function(object) {
                    if(!is.na(object@.w))
                      stop("[Validation: power] .w is an internal variable that should not be set by user")
                    }, contains = "Elementwise")
setMethod("initialize", "power", function(.Object, ..., x, p, max_denom = 1024, .w = NA_real_) {
  # TODO: Fill in p and w accordingly
  .Object@x <- x
  .Object@max_denom <- max_denom
  callNextMethod(.Object, ..., .args = list(x))
})

setMethod("sign_from_args", "power", function(object) {
  if(object@p == 1)
    object@.args[1]@dcp_attr@sign
  else
    Sign(sign = SIGN_POSITIVE_KEY)
})

setMethod("func_curvature", "power", function(object) {
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

setMethod("monotonicity", "power", function(object) {
  if(object@p ==0)
    list(INCREASING)
  else if(object@p == 1)
    list(INCREASING)
  else if(object@p < 0)
    list(DECREASING)
  else if(object@p > 0 && object@p < 1)
    list(INCREASING)
  else if(object@p > 1) {
    if(is_power2(object@p))
      list(SIGNED)
    else
      list(INCREASING)
  }
  else
    stop("Unknown monotonicity for power p = ", object@p)
})

setMethod("get_data", "power", function(object) { list(object@p, object@w) })
