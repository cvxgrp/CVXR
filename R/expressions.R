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
setMethod("curvature", "Expression", function(object) {
  if(is_constant(object))
    curvature_str <- Curvature.CONSTANT
  else if(is_affine(object))
    curvature_str <- Curvature.AFFINE
  else if(is_convex(object))
    curvature_str <- Curvature.CONVEX
  else if(is_concave(object))
    curvature_str <- Curvature.CONCAVE
  else
    curvature_str <- Curvature.UNKNOWN
  curvature_str
})
setMethod("is_constant", "Expression", function(object) { length(variables(object)) == 0 || is_zero(object) })
setMethod("is_affine", "Expression", function(object) { is_constant(object) || (is_convex(object) && is_concave(object)) })
setMethod("is_convex", "Expression", function(object) { stop("Unimplemented") })
setMethod("is_concave", "Expression", function(object) { stop("Unimplemented") })
setMethod("is_dcp", "Expression", function(object) { is_convex(object) || is_concave(object) })
setMethod("is_quadratic", "Expression", function(object) { FALSE })

# Sign properties
setMethod("sign", "Expression", function(object) {
  if(is_zero(object))
    sign_str <- Sign.ZERO
  else if(is_positive(object))
    sign_str <- Sign.POSITIVE
  else if(is_negative(object))
    sign_str <- Sign.NEGATIVE
  else
    sign_str <- Sign.UNKNOWN
  sign_str
})
setMethod("is_zero", "Expression", function(object) { is_positive(object) && is_negative(object) })
setMethod("is_positive", "Expression", function(object) { stop("Unimplemented") })
setMethod("is_negative", "Expression", function(object) { stop("Unimplemented") })
setMethod("size", "Expression", function(object) { stop("Unimplemented") })
setMethod("is_scalar", "Expression", function(object) { all(size(self) == c(1,1)) })
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
setMethod("^", signature(e1 = "Expression", e2 = "numeric"), function(e1, e2) { Power(x = e1, p = e2) })

# Matrix operators
t.Expression <- function(x) { if(is_scalar(x)) x else Transpose(.args = list(x)) }   # Need S3 method dispatch as well
setMethod("t", signature(x = "Expression"), function(x) { if(is_scalar(x)) x else Transpose(.args = list(x)) })
# TODO: Overload the [ operator for slicing rows/columns from an expression
setMethod("%*%", signature(x = "Expression", y = "Expression"), function(x, y) {
  # Multiplying by a constant on the right is handled differently
  # from multiplying by a constant on the left
  if(is_constant(x)) {
    if(size(x)[1] == size(y)[1] && size(x)[2] != size(y)[1] && is(x, "Constant") && is_1d_array(x))
      x <- t(x)
    return(MulExpression(lh_exp = x, rh_exp = y))
  } else if(is_constant(y)) {
    # Having the constant on the left is more efficient
    if(is_scalar(x) || is_scalar(y))
      return(MulExpression(lh_exp = y, rh_exp = x))
    else
      return(RMulExpression(lh_exp = x, rh_exp = y))
  # When both expressions are not constant, allow affine * affine, but raise DCPError otherwise
  # Cannot multiply two non-constant expressions
  } else if(is_affine(x) && is_affine(y)) {
    warn("Forming a non-convex expression (affine) * (affine)")
    return(AffineProd(x = x, y = y))
  } else
    stop("Cannot multiply ", curvature(x), "and", curvature(y))
})
setMethod("%*%", signature(x = "Expression", y = "ConstVal"), function(x, y) { x %*% as.Constant(y) })
setMethod("%*%", signature(x = "ConstVal", y = "Expression"), function(x, y) { as.Constant(x) %*% y })

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

# Positive definite inequalities
setMethod("%>>%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e1, e2) })
setMethod("%>>%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %>>% as.Constant(e2) })
setMethod("%>>%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %>>% e2 })
setMethod("%<<%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e2, e1) })
setMethod("%<<%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %<<% as.Constant(e2) })
setMethod("%<<%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %<<% e2 })

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
setMethod("is_quadratic", "Leaf", function(object) { TRUE })
setMethod("domain", "Leaf", function(object) { list() })   # Default is full domain

# TODO: Check this is doing the correct thing
setMethod("validate_val", "Leaf", function(object, val) {
  if(!is.na(val)) {
    # Convert val to the proper matrix type
    val <- as.matrix(val)
    dims <- dim(val)
    if(dims != size(object))
      stop("Invalid dimensions for value")
    
    # All signs are valid if sign is unknown
    # Otherwise value sign must match declared sign
    pos_val <- min(val) >= -1e-5
    neg_val <- max(val) <= 1e-5
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
