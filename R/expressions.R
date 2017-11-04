#'
#' The Expression class.
#'
#' This class represents a mathematical expression.
#'
#' @name Expression
#' @rdname Expression-class
#' @export
Expression <- setClass("Expression", contains = "Canonical")

setOldClass("data.frame")
setOldClass("matrix")
setOldClass("vector")
setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))
setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "vector", "numeric", "dMatrix"))
setClassUnion("ConstValORExpr", c("ConstVal", "Expression"))
setClassUnion("ListORExpr", c("list", "Expression"))
setClassUnion("ConstValListORExpr", c("ConstVal", "list", "Expression"))
setClassUnion("NumORgmp", c("numeric", "bigq", "bigz"))

# Helper function since syntax is different for LinOp (list) vs. Expression object
setMethod("size", "ListORExpr", function(object) {
  if(is.list(object))
    object$size
  else
    size(object)
})

# Casts the second argument of a binary operator as an Expression
.cast_other <- function(binary_op) {
  cast_op <- function(object, other) {
    other <- as.Constant(other)
    binary_op(object, other)
  }
  cast_op
}

#' @rdname value-methods
#' @describeIn Expression The value of the expression.
setMethod("value", "Expression", function(object) { stop("Unimplemented") })

#' @rdname grad
#' @describeIn Expression The (sub/super)-gradient of the expression with respect to each variable.
setMethod("grad", "Expression", function(object) { stop("Unimplemented") })

#' @rdname domain
#' @describeIn Expression A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "Expression", function(object) { stop("Unimplemented") })

setMethod("show", "Expression", function(object) {
  cat("Expression(", curvature(object), ", ", sign(object), ", ", size(object), ")", sep = "")
})

# A string with information about the expression such as curvature, sign, and size.
setMethod("as.character", "Expression", function(x) {
  paste("Expression(", curvature(x), ", ", sign(x), ", ", size(x), ")", sep = "")
})

#' @rdname name
#' @describeIn Expression The string representation of the expression.
setMethod("name", "Expression", function(object) { stop("Unimplemented") })

#' @rdname curvature
#' @describeIn Expression The curvature of the expression.
setMethod("curvature", "Expression", function(object) {
  if(is_constant(object))
    curvature_str <- CONSTANT
  else if(is_affine(object))
    curvature_str <- AFFINE
  else if(is_convex(object))
    curvature_str <- CONVEX
  else if(is_concave(object))
    curvature_str <- CONCAVE
  else
    curvature_str <- UNKNOWN
  curvature_str
})

#' @rdname curvature-methods
#' @describeIn Expression The expression is constant if it contains no variables or is identically zero.
setMethod("is_constant", "Expression", function(object) { length(variables(object)) == 0 || is_zero(object) })

#' @rdname curvature-methods
#' @describeIn Expression The expression is affine if it is constant or both convex and concave.
setMethod("is_affine", "Expression", function(object) { is_constant(object) || (is_convex(object) && is_concave(object)) })

#' @rdname curvature-methods
#' @describeIn Expression A logical value indicating whether the expression is convex.
setMethod("is_convex", "Expression", function(object) { stop("Unimplemented") })

#' @rdname curvature-methods
#' @describeIn Expression A logical value indicating whether the expression is concave.
setMethod("is_concave", "Expression", function(object) { stop("Unimplemented") })

#' @rdname is_dcp
#' @describeIn Expression The expression is DCP if it is convex or concave.
setMethod("is_dcp", "Expression", function(object) { is_convex(object) || is_concave(object) })

#' @rdname curvature-methods
#' @describeIn Expression A logical value indicating whether the expression is quadratic.
setMethod("is_quadratic", "Expression", function(object) { FALSE })

#' @rdname curvature-methods
#' @describeIn Expression A logical value indicating whether the expression is piecewise linear.
setMethod("is_pwl", "Expression", function(object) { FALSE })

#'
#' Sign of Expression
#' 
#' The sign of an expression.
#' 
#' @param x An \linkS4class{Expression} object.
#' @return A string indicating the sign of the expression, either "ZERO", "POSITIVE", "NEGATIVE", or "UNKNOWN".
#' @docType methods
#' @rdname sign
#' @export
setMethod("sign", "Expression", function(x) {
  if(is_zero(x))
    sign_str <- ZERO
  else if(is_positive(x))
    sign_str <- POSITIVE
  else if(is_negative(x))
    sign_str <- NEGATIVE
  else
    sign_str <- UNKNOWN
  sign_str
})

#' @rdname sign-methods
#' @describeIn Expression The expression is zero if it is both positive and negative.
setMethod("is_zero", "Expression", function(object) { is_positive(object) && is_negative(object) })

#' @rdname sign-methods
#' @describeIn Expression A logical value indicating whether the expression is positive.
setMethod("is_positive", "Expression", function(object) { stop("Unimplemented") })

#' @rdname sign-methods
#' @describeIn Expression A logical value indicating whether the expression is negative.
setMethod("is_negative", "Expression", function(object) { stop("Unimplemented") })

#' @rdname size
#' @describeIn Expression The \code{c(row, col)} dimensions of the expression.
setMethod("size", "Expression", function(object) { stop("Unimplemented") })

#' @rdname size-methods
#' @describeIn Expression The expression is scalar if rows = cols = 1.
setMethod("is_scalar", "Expression", function(object) { all(size(object) == c(1,1)) })

#' @rdname size-methods
#' @describeIn Expression The expression is a vector if min(rows, cols) = 1.
setMethod("is_vector", "Expression", function(object) { min(size(object)) == 1 })

#' @rdname size-methods
#' @describeIn Expression The expression is a matrix if rows > 1 and cols > 1.
setMethod("is_matrix", "Expression", function(object) { size(object)[1] > 1 && size(object)[2] > 1 })

#' @describeIn Expression Number of rows in the expression.
#' @export
setMethod("nrow", "Expression", function(x) { size(x)[1] })

#' @describeIn Expression Number of columns in the expression.
#' @export
setMethod("ncol", "Expression", function(x) { size(x)[2] })

# Slice operators
#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "missing", j = "missing", drop = "ANY"), function(x, i, j, ..., drop) { x })

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "index", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  if(is_vector(x) && size(x)[1] < size(x)[2])
    Index.get_special_slice(x, NULL, i)   # If only first index given, apply it along longer dimension of vector
  else
    Index.get_special_slice(x, i, NULL)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "missing", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, NULL, j)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "index", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "matrix", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "index", j = "matrix", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "matrix", j = "matrix", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "matrix", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  # This follows conventions in Matrix package, but differs from base handling of matrices
  Index.get_special_slice(x, i, NULL)
})

#' @rdname Index-class
setMethod("[", signature(x = "Expression", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  stop("Invalid or unimplemented Expression slice operation")
})

# Arithmetic operators
#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { e1 })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { NegExpression(expr = e1) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e2 + e1 })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 + NegExpression(expr = e2) })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 + (-e2) })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e1 + NegExpression(expr = e2) })

#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(is_constant(e1))
    MulElemwise(lh_const = e1, rh_exp = e2)
  else if(is_constant(e2))
    MulElemwise(lh_const = e2, rh_exp = e1)
  else if(all(size(e1) == c(1,1)) && all(size(e2) == c(1,1)) && is_affine(e1) && is_affine(e2)) {
    warning("Forming a non-convex expression (affine) * (affine)")
    AffineProd(x = e1, y = e2)
  } else
    stop("Cannot multiply elementwise ", curvature(e1), " and ", curvature(e2))
})

#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { as.Constant(e2) * e1 })

#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) * e2 })

#' @rdname DivExpression-class
setMethod("/", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(is_constant(e2) && is_scalar(e2))
    DivExpression(lh_exp = e1, rh_exp = e2)
  else
    stop("Can only divide by a scalar constant")
})

#' @rdname DivExpression-class
setMethod("/", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 / as.Constant(e2) })

#' @rdname DivExpression-class
setMethod("/", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) / e2 })

#' @docType methods
#' @rdname power
setMethod("^", signature(e1 = "Expression", e2 = "numeric"), function(e1, e2) {
  if(e2 == 2)
    Square(x = e1)
  else if(e2 == 0.5)
    Sqrt(x = e1)
  else
    Power(x = e1, p = e2)
})

# Matrix operators
#'
#' Matrix Transpose
#'
#' The transpose of a matrix.
#'
#' @param x An \linkS4class{Expression} representing a matrix.
#' @return An \linkS4class{Expression} representing the transposed matrix.
#' @docType methods
#' @rdname transpose
#' @aliases t
#' @export
t.Expression <- function(x) { if(is_scalar(x)) x else Transpose(args = list(x)) }   # Need S3 method dispatch as well

#' @docType methods
#' @rdname transpose
#' @export
setMethod("t", signature(x = "Expression"), function(x) { if(is_scalar(x)) x else Transpose(args = list(x)) })

#' @rdname MulExpression-class
setMethod("%*%", signature(x = "Expression", y = "Expression"), function(x, y) {
  # if(is_scalar(x) || is_scalar(y))
  #  stop("Scalar operands are not allowed, use '*' instead")
  
  # Multiplying by a constant on the right is handled differently
  # from multiplying by a constant on the left
  if(is_constant(x)) {
    if(size(x)[1] == size(y)[1] && size(x)[2] != size(y)[1] && is(x, "Constant") && x@is_1D_array)
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
    warning("Forming a non-convex expression (affine) * (affine)")
    return(AffineProd(x = x, y = y))
  } else
    stop("Cannot multiply ", curvature(x), " and ", curvature(y))
})

#' @rdname MulExpression-class
setMethod("%*%", signature(x = "Expression", y = "ConstVal"), function(x, y) { x %*% as.Constant(y) })

#' @rdname MulExpression-class
setMethod("%*%", signature(x = "ConstVal", y = "Expression"), function(x, y) { as.Constant(x) %*% y })

# Comparison operators
#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { EqConstraint(lh_exp = e1, rh_exp = e2) })

#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 == as.Constant(e2) })

#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) == e2 })

#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { LeqConstraint(lh_exp = e1, rh_exp = e2) })

#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 <= as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) <= e2 })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 <= e2 })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 < as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) < e2 })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e2 <= e1 })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 >= as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) >= e2 })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 >= e2 })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 > as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) > e2 })

# Positive definite inequalities
#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%>>%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e1, e2) })

#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%>>%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %>>% as.Constant(e2) })

#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%>>%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %>>% e2 })

#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%<<%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e2, e1) })

#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%<<%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %<<% as.Constant(e2) })

#' @docType methods
#' @rdname PSDConstraint
#' @export
setMethod("%<<%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %<<% e2 })

#'
#' The Leaf class.
#'
#' This class represents a leaf node, i.e. a Variable, Constant, or Parameter.
#'
#' @slot args A list containing the arguments.
#' @name Leaf
#' @rdname Leaf-class
Leaf <- setClass("Leaf", representation(args = "list"), prototype(args = list()), contains = "Expression")

#' @rdname expression-parts
#' @describeIn Leaf List of \linkS4class{Variable} objects in the leaf node.
setMethod("variables", "Leaf", function(object) { list() })

#' @rdname expression-parts
#' @describeIn Leaf List of \linkS4class{Parameter} objects in the leaf node.
setMethod("parameters", "Leaf", function(object) { list() })

#' @rdname expression-parts
#' @describeIn Leaf List of \linkS4class{Constant} objects in the leaf node.
setMethod("constants", "Leaf", function(object) { list() })

#' @rdname curvature-methods
#' @describeIn Leaf A logical value indicating whether the leaf node is convex.
setMethod("is_convex", "Leaf", function(object) { TRUE })

#' @rdname curvature-methods
#' @describeIn Leaf A logical value indicating whether the leaf node is concave.
setMethod("is_concave", "Leaf", function(object) { TRUE })

#' @rdname curvature-methods
#' @describeIn Leaf A logical value indicating whether the leaf node is quadratic.
setMethod("is_quadratic", "Leaf", function(object) { TRUE })

#' @rdname curvature-methods
#' @describeIn Leaf A logical value indicating whether the leaf node is piecewise linear.
setMethod("is_pwl", "Leaf", function(object) { TRUE })

#' @rdname domain
#' @describeIn A list of constraints describing the closure of the region where the leaf node is finite. Default is the full domain.
setMethod("domain", "Leaf", function(object) { list() })   # Default is full domain

setMethod("validate_val", "Leaf", function(object, val) {
  if(length(val) > 1 || !(length(val) == 1 && is.na(val))) {
    # Convert val to the proper matrix type
    if(!(is.null(dim(val)) && length(val) == 1))
      val <- as.matrix(val)
    size <- intf_size(val)
    if(any(size != size(object)))
      stop("Invalid dimensions (", size[1], ", ", size[2], ") for ", class(object), " value")
    
    # All signs are valid if sign is unknown
    # Otherwise value sign must match declared sign
    sign <- intf_sign(val)
    pos_val <- sign[1]
    neg_val <- sign[2]
    if(is_positive(object) && !pos_val || is_negative(object) && !neg_val)
      stop("Invalid sign for ", class(object), " value")
    # Round to correct sign
    else if(is_positive(object))
      val <- pmax(val, 0)
    else if(is_negative(object))
      val <- pmin(val, 0)
  }
  return(val)
})
