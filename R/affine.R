AffAtom <- setClass("AffAtom", contains = "Atom")
setMethod("func_curvature", "AffAtom", function(object) { Curvature(CURV_AFFINE_KEY) })
setMethod("sign_from_args", "AffAtom", function(object) { 
  arg_signs = lapply(get_slots(object), function(x) { x@dcp_attr@sign })
  Reduce("+", arg_signs)
})
setMethod("monotonicity", "AffAtom", function(object) {
  rep(INCREASING, length(object@.args))
})

AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")
setMethod("init_dcp_attr", "AddExpression", function(object) {
  arg_dcp <- lapply(object@.args, function(x) { x@dcp_attr })
  Reduce("+", arg_dcp)
})
setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., .args = arg_groups)   # Casts R values to Constant objects
  .Object@.args <- lapply(arg_groups, function(group) { if(is(group,"AddExpression")) group@.args else group })
  .Object@.args <- flatten_list(.Object@.args)   # Need to flatten list of expressions
  return(.Object)
})

UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression", op_name = "character"), contains = "AffAtom")
setMethod("initialize", "UnaryOperator", function(.Object, ..., expr, op_name) {
  .Object@expr = expr
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})
setMethod("init_dcp_attr", "UnaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr)
})

NegExpression <- setClass("NegExpression", contains = "UnaryOperator")
setMethod("initialize", "NegExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "-")
})
NegExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(neg_expr(arg_objs[1]), list())
}

BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "Expression", rh_exp = "Expression", op_name = "character"), contains = "AffAtom")
setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp, op_name) {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@lh_exp, .Object@rh_exp))
})
setMethod("init_dcp_attr", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr, object@.args[[2]]@dcp_attr)
})
setMethod("validate_args", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr@shape, object@.args[[2]]@dcp_attr@shape)
})

MulExpression <- setClass("MulExpression", contains = "BinaryOperator")
setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "*")
})
MulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[2] != 1 && all(size(arg_objs[2]) == c(1,1))) {
    arg <- promote(arg_objs[2], list(size(2), 1))
    arg_objs[2] <- diag_vec(arg)
  }
  list(mul_expr(arg_objs[1], arg_objs[2], size), list())
}

RMulExpression <- setClass("RMulExpression", contains = "MulExpression")
RMulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[1] != 1 && all(size(arg_objs[1]) == c(1,1))) {
    arg <- promote(arg_objs[1], list(size[1], 1))
    arg_objs[1] <- diag_vec(arg)
  }
  list(rmul_expr(arg_objs[1], arg_objs[2], size), list())
}

DivExpression <- setClass("DivExpression", contains = "BinaryOperator")
setMethod("initialize", "DivExpression", function(.Object, ...) {
  .Object@op_name = "/"
  callNextMethod(.Object, ..., op_name = "/")
})
DivExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(div_expr(arg_objs[1], arg_objs[2]), list())
}

Transpose <- setClass("Transpose", contains = "AffAtom")
setMethod("shape_from_args", "Transpose", function(object) {
  obj_size = size(object@.args[[1]])
  Shape(rows = obj_size[2], cols = obj_size[1])
})
