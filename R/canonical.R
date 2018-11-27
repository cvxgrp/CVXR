#'
#' The Canonical class.
#'
#' This virtual class represents a canonical expression.
#'
#' @rdname Canonical-class
setClass("Canonical", contains = "VIRTUAL")

#' @param object A \linkS4class{Canonical} object.
#' @describeIn Canonical The expression associated with the input.
setMethod("expr", "Canonical", function(object) {
  if(length(object@args) != 1)
    stop("'expr' is ambiguous, there should only be one argument.")
  return(object@args[[1]])
})

#' @describeIn Canonical The graph implementation of the input.
setMethod("canonicalize", "Canonical", function(object) { stop("Unimplemented") })
setMethod("canonical_form", "Canonical", function(object) { canonicalize(object) })

#' @describeIn Canonical List of \linkS4class{Variable} objects in the expression.
setMethod("variables", "Canonical", function(object) {
  unique(flatten_list(lapply(object@args, function(arg) { variables(arg) })))
})

#' @describeIn Canonical List of \linkS4class{Parameter} objects in the expression.
setMethod("parameters", "Canonical", function(object) {
  unique(flatten_list(lapply(object@args, function(arg) { parameters(arg) })))
})

#' @describeIn Canonical List of \linkS4class{Constant} objects in the expression.
setMethod("constants", "Canonical", function(object) {
  const_list <- flatten_list(lapply(object@args, function(arg) { constants(arg) }))
  const_id <- sapply(const_list, function(constant) { id(constant) })
  const_list[!duplicated(const_id)]
})

#' @describeIn Canonical List of \linkS4class{Atom} objects in the expression.
setMethod("atoms", "Canonical", function(object) {
  unique(flatten_list(lapply(object@args, function(arg) { atoms(arg) })))
})

#' @describeIn Canonical Information needed to reconstruct the expression aside from its arguments.
setMethod("get_data", "Canonical", function(object) { list() })

