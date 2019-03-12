#'
#' The Canonical class.
#'
#' This virtual class represents a canonical expression.
#'
#' @rdname Canonical-class
setClass("Canonical", representation(args = "list"), prototype(args = list()), contains = "VIRTUAL")

#' @param object A \linkS4class{Canonical} object.
#' @describeIn Canonical The expression associated with the input.
setMethod("expr", "Canonical", function(object) {
  if(length(object@args) != 1)
    stop("'expr' is ambiguous, there should only be one argument.")
  return(object@args[[1]])
})

#' @describeIn Canonical The graph implementation of the input.
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
  # const_id <- sapply(const_list, function(constant) { id(constant) })
  # const_list[!duplicated(const_id)]
  const_list[!duplicated(const_list)]
})

#' @describeIn Canonical List of \linkS4class{Atom} objects in the expression.
setMethod("atoms", "Canonical", function(object) {
  unique(flatten_list(lapply(object@args, function(arg) { atoms(arg) })))
})

setMethod("tree_copy", "Canonical", function(object, id_objects = list()) {
  new_args <- list()
  for(arg in object@args) {
    if(is.list(arg)) {
      arg_list <- lapply(arg, function(elem) { tree_copy(elem, id_objects) })
      new_args <- c(new_args, arg_list)
    } else
      new_args <- c(new_args, tree_copy(arg, id_objects))
  }
  return(copy(object, args = new_args, id_objects = id_objects))
})

setMethod("copy", "Canonical", function(object, args = NULL, id_objects = list()) {
  if("id" %in% names(attributes(object)) && as.character(object@id) %in% names(id_objects))
    return(id_objects[[as.character(object@id)]])
  if(is.null(args))
    args <- object@args
  data <- get_data(object)
  if(!is.null(data) && length(data) != 0)
    return(do.call(class(object), c(args, data)))
  else
    return(do.call(class(object), args))
})

#' @describeIn Canonical Information needed to reconstruct the expression aside from its arguments.
setMethod("get_data", "Canonical", function(object) { list() })

