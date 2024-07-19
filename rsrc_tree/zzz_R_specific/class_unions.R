setOldClass("data.frame")
setOldClass("matrix")
setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))

setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "numeric", "complex", "dMatrix", "bigq", "bigz"))
setClassUnion("ConstValORExpr", c("ConstVal", "Expression"))
setClassUnion("ConstValORNULL", c("ConstVal", "NULL"))
setClassUnion("ConstValListORExpr", c("ConstVal", "list", "Expression"))
setClassUnion("ListORExpr", c("list", "Expression"))
setClassUnion("NumORgmp", c("numeric", "bigq", "bigz"))
setClassUnion("NumORNULL", c("numeric", "NULL"))
setClassUnion("NumORLogical", c("logical", "numeric"))
#setClassUnion("ReducedMatORNULL", c("ReducedMat", "NULL"))
setClassUnion("S4ORNULL", c("S4", "NULL"))

# Helper function since syntax is different for LinOp (list) vs. Expression object
#' @rdname size
setMethod("size", "ListORExpr", function(object) {
  if(is.list(object))
    object$size
  else
    size(object)
})

# Helper function so we can flatten both Expression objects and regular matrices into a single column vector.
setMethod("flatten", "numeric", function(object) { matrix(object, ncol = 1) })

