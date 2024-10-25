
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

