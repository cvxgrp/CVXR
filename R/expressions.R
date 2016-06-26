#'
#' The Expression class.
#'
#' This class represents an expression in CVXR.
#'
#' @slot dcp_attr A \code{DCPAttr} object specifying the DCP attributes of the expression.
#' @aliases Expression
#' @export
Expression <- setClass("Expression", representation(dcp_attr = "DCPAttr"), prototype(dcp_attr = new("DCPAttr")))

setOldClass("data.frame")
setOldClass("matrix")
setOldClass("vector")
setClassUnion("ConstVal", c("data.frame", "matrix", "vector", "numeric"))
setClassUnion("ConstValORExpr", c("ConstVal", "Expression"))

setMethod("show", "Expression", function(object) {
  cat("Expression(", as.character(object@dcp_attr@curvature), ", ", as.character(object@dcp_attr@sign), ", ", as.character(object@dcp_attr@shape), ")", sep = "")
})

setMethod("as.character", "Expression", function(x) {
  paste("Expression(", as.character(x@dcp_attr@curvature), ", ", as.character(x@dcp_attr@sign), ", ", as.character(x@dcp_attr@shape), ")", sep = "")
})

setMethod("canonical_form", "Expression", function(object) { canonicalize(object) })
setMethod("get_data", "Expression", function(object) { list() })

# Curvature properties
setMethod("curvature", "Expression", function(object) { object@dcp_attr@curvature })
setMethod("is_constant", "Expression", function(object) { is_constant(object@dcp_attr@curvature) })
setMethod("is_affine", "Expression", function(object) { is_affine(object@dcp_attr@curvature) })
setMethod("is_convex", "Expression", function(object) { is_convex(object@dcp_attr@curvature) })
setMethod("is_concave", "Expression", function(object) { is_concave(object@dcp_attr@curvature) })
setMethod("is_dcp", "Expression", function(object) { is_dcp(object@dcp_attr@curvature) })

# Sign properties
setMethod("sign", "Expression", function(x) { x@dcp_attr@sign })
setMethod("is_zero", "Expression", function(object) { is_zero(object@dcp_attr@sign) })
setMethod("is_positive", "Expression", function(object) { is_positive(object@dcp_attr@sign) })
setMethod("is_negative", "Expression", function(object) { is_negative(object@dcp_attr@sign) })

# Shape properties
setMethod("size", "Expression", function(object) { size(object@dcp_attr@shape) })
setMethod("is_scalar", "Expression", function(object) { all(size(object) == c(1,1)) })
setMethod("is_vector", "Expression", function(object) { min(size(object)) == 1 })
setMethod("is_matrix", "Expression", function(object) { size(object)[1] > 1 && size(object)[2] > 1 })

# Arithmetic operators
setMethod("+", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { e1 })
setMethod("-", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { NegExpression(expr = e1) })
setMethod("+", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })
setMethod("+", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })
setMethod("+", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e2 + e1 })
setMethod("-", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 + -e2 })
setMethod("-", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 + -e2 })
setMethod("-", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e2 - e1 })
setMethod("*", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(!is_constant(e1) && !is_constant(e2))
    stop("Cannot multiply two non-constants")
  else if(is_constant(e1)) {
    if(size(e1)[1] == size(e2)[1] && size(e1)[2] != size(e2)[1] && is(e1, "Constant"))
      MulExpression(lh_exp = t(e1), rh_exp = e2)
    else
      MulExpression(lh_exp = e1, rh_exp = e2)
  } else if(is_scalar(e1) || is_scalar(e2))
    MulExpression(lh_exp = e2, rh_exp = e1)
  else
    RMulExpression(lh_exp = e1, rh_exp = e2)
})
setMethod("*", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 * as.Constant(e2) })
setMethod("*", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) * e2 })
setMethod("/", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(is_constant(e2) && is_scalar(e2))
    DivExpression(lh_exp = e1, rh_exp = e2)
  else
    stop("Can only divide by a scalar constant")
})
setMethod("/", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 / as.Constant(e2) })
setMethod("/", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) / e2 })

# Comparison operators
setMethod("==", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { EqConstraint(lh_exp = e1, rh_exp = e2) })
setMethod("==", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 == as.Constant(e2) })
setMethod("==", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) == e2 })
setMethod("<=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { LeqConstraint(lh_exp = e1, rh_exp = e2) })
setMethod("<=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 <= as.Constant(e2) })
setMethod("<=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) <= e2 })
setMethod("<",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 <= e2 })
setMethod("<",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 < as.Constant(e2) })
setMethod("<",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) < e2 })
setMethod(">=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e2 <= e1 })
setMethod(">=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 >= as.Constant(e2) })
setMethod(">=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) >= e2 })
setMethod(">",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 >= e2 })
setMethod(">",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 > as.Constant(e2) })
setMethod(">",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) > e2 })

t.Expression <- function(x) { if(is_scalar(x)) x else Transpose(.args = list(x)) }   # Need S3 method dispatch as well
setMethod("t", signature(x = "Expression"), function(x) { if(is_scalar(x)) x else Transpose(.args = list(x)) })
setMethod("^", signature(e1 = "Expression", e2 = "numeric"), function(e1, e2) { Power(x = e1, p = e2) })
# TODO: Should I overload matrix multiplication %*% operator to point to regular multiplication *?
# TODO: Overload the [ operator for slicing rows/columns from an expression

#'
#' The Leaf class.
#'
#' This class represents a leaf node, i.e. a Variable, Constant, or Parameter.
#'
Leaf <- setClass("Leaf", contains = "Expression")

setMethod("variables", "Leaf", function(object) { list() })
setMethod("parameters", "Leaf", function(object) { list() })

