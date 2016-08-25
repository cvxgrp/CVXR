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
setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))
setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "vector", "numeric"))
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
setMethod("*", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { MulElemwise(lh_const = e1, rh_exp = e2) })
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
setMethod("%*%", signature(x = "Expression", y = "Expression"), function(x, y) {
  if(!is_constant(x) && !is_constant(y))
    stop("Cannot multiply two non-constants")
  else if(is_constant(x)) {
    if(size(x)[1] == size(y)[1] && size(x)[2] != size(y)[1] && is(x, "Constant"))
      MulExpression(lh_exp = t(x), rh_exp = y)
    else
      MulExpression(lh_exp = x, rh_exp = y)
  } else if(is_scalar(x) || is_scalar(y))
    MulExpression(lh_exp = y, rh_exp = x)
  else
    RMulExpression(lh_exp = x, rh_exp = y)
})
setMethod("%*%", signature(x = "Expression", y = "ConstVal"), function(x, y) { x %*% as.Constant(y) })
setMethod("%*%", signature(x = "ConstVal", y = "Expression"), function(x, y) { as.Constant(x) %*% y })
# TODO: Overload the [ operator for slicing rows/columns from an expression

#'
#' The Leaf class.
#'
#' This class represents a leaf node, i.e. a Variable, Constant, or Parameter.
#'
Leaf <- setClass("Leaf", contains = "Expression")

setMethod("variables", "Leaf", function(object) { list() })
setMethod("parameters", "Leaf", function(object) { list() })
setMethod("constants", "Leaf", function(object) { list() })
setMethod("is_convex", "Leaf", function(object) { TRUE })
setMethod("is_concave", "Leaf", function(object) { TRUE })
setMethod("domain", "Leaf", function(object) { list() })
setMethod("validate_value", "Leaf", function(object, val) { 
  if(!is.na(val)) {
    # Convex val to the proper matrix type
    val <- as.matrix(val)
    dims <- dim(val)
    if(dims != size(object))
      stop("Invalid dimensions for value")
    
    # All signs are valid if sign is unknown
    # Otherwise value sign must match declared sign
    pos_val <- min(val) >= 0
    neg_val <- max(val) <= 0
    if(is_positive(object) && !pos_val || is_negative(object) && !neg_val)
      stop("Invalid sign for value")
    # Round to correct sign
    else if(is_positive(object))
      val <- max(val, 0)
    else if(is_negative(object))
      val <- min(val, 0)
  }
  return(val)
})
