#' @importFrom utils assignInMyNamespace
#' @importFrom R.utils intToBin
#' @importFrom Rmpfr mpfr getPrec
#' @importFrom gmp as.bigq as.bigz is.bigq is.bigz is.whole numerator denominator asNumeric
#' @importFrom bit64 as.integer64 as.bitstring
#' @importClassesFrom gmp bigq bigz

##
## Global variable for package CVXR
## idCounter, numpy handle, scipy.sparse handle
##
.CVXR.options <- list(idCounter = 0L, np = NULL, sp = NULL, mosekglue = NULL)

#'
#' Set ID Counter
#'
#' Set the CVXR variable/constraint identification number counter.
#'
#' @param value The value to assign as ID.
#' @return the changed value of the package global \code{.CVXR.options}.
#' @export
#' @examples
#' \dontrun{
#'   setIdCounter(value = 0L)
#' }
setIdCounter <- function(value = 0L) {
    .CVXR.options$idCounter <- value
    assignInMyNamespace(".CVXR.options", .CVXR.options)
    .CVXR.options
}

#'
#' Reset Options
#'
#' Reset the global package variable \code{.CVXR.options}.
#'
#' @return The default value of CVXR package global \code{.CVXR.options}.
#' @export
#' @examples
#' \dontrun{
#'   resetOptions()
#' }
resetOptions <- function() {
    assignInMyNamespace(".CVXR.options", list(idCounter = 0L, np = NULL, sp = NULL, mosekglue = NULL))
    .CVXR.options
}

#'
#' Get ID
#'
#' Get the next identifier value.
#'
#' @return A new unique integer identifier.
#' @export
#' @examples
#' \dontrun{
#'    get_id()
#' }
get_id <- function() {
    id <- .CVXR.options$idCounter <- .CVXR.options$idCounter + 1L
    assignInMyNamespace(".CVXR.options", .CVXR.options)
    id
}

#'
#' Get scipy handle
#'
#' Get the scipy handle or fail if not available
#'
#' @return the scipy handle
#' @export
#' @examples
#' \dontrun{
#'    get_sp
#' }
get_sp <- function() {
    sp <- .CVXR.options$sp
    if (is.null(sp)) {
        stop("Scipy not available")
    }
    sp
}

#'
#' Get numpy handle
#'
#' Get the numpy handle or fail if not available
#'
#' @return the numpy handle
#' @export
#' @examples
#' \dontrun{
#'    get_np
#' }
get_np <- function() {
    np <- .CVXR.options$np
    if (is.null(np)) {
        stop("Numpy not available")
    }
    np
}

#'
#' Get our mosekglue handle
#'
#' Get the mosekglue handle or fail if not available
#'
#' @return the mosekglue handle
#' @export
#' @examples
#' \dontrun{
#'    get_mosekglue
#' }
get_mosekglue <- function() {
    mosekglue <- .CVXR.options$mosekglue
    if (is.null(mosekglue)) {
        stop("CVXR python mosekglue not available")
    }
    mosekglue
}

#'
#' Get our gurobiglue handle
#'
#' Get the gurobiglue handle or fail if not available
#'
#' @return the gurobiglue handle
#' @export
#' @examples
#' \dontrun{
#'    get_gurobiglue
#' }
get_gurobiglue <- function() {
    gurobiglue <- .CVXR.options$gurobiglue
    if (is.null(gurobiglue)) {
        stop("CVXR python gurobiglue not available")
    }
    gurobiglue
}

#'
#' The Rdict class.
#'
#' A simple, internal dictionary composed of a list of keys and a list of values. These keys/values can be any type, including nested lists, S4 objects, etc.
#' Incredibly inefficient hack, but necessary for the geometric mean atom, since it requires mixed numeric/gmp objects.
#'
#' @slot keys A list of keys.
#' @slot values A list of values corresponding to the keys.
#' @name Rdict-class
#' @aliases Rdict
#' @rdname Rdict-class
setClass("Rdict", representation(keys = "list", values = "list"), prototype(keys = list(), values = list()),
         validity = function(object) {
           if(length(object@keys) != length(object@values))
             return("Number of keys must match number of values")
           if(!all(unique(object@keys) != object@keys))
             return("Keys must be unique")
           return(TRUE)
         })

#' @param keys A list of keys.
#' @param values A list of values corresponding to the keys.
#' @rdname Rdict-class
Rdict <- function(keys = list(), values = list()) {
  new("Rdict", keys = keys, values = values)
}

#' @param x,set A \linkS4class{Rdict} object.
#' @param name Either "keys" for a list of keys, "values" for a list of values, or "items" for a list of lists where each nested list is a (key, value) pair.
#' @rdname Rdict-class
setMethod("$", signature(x = "Rdict"), function(x, name) {
  if(name == "items") {
    items <- rep(list(list()), length(x))
    for(i in 1:length(x)) {
      tmp <- list(key = x@keys[[i]], value = x@values[[i]])
      items[[i]] <- tmp
    }
    return(items)
  } else
    slot(x, name)
})

#' @rdname Rdict-class
setMethod("length", signature(x = "Rdict"), function(x) { length(x@keys) })

#' @param el The element to search the dictionary of values for.
#' @rdname Rdict-class
setMethod("is.element", signature(el = "ANY", set = "Rdict"), function(el, set) {
  for(k in set@keys) {
    if(identical(k, el))
      return(TRUE)
  }
  return(FALSE)
})

#' @param i A key into the dictionary.
#' @param j,drop,... Unused arguments.
#' @rdname Rdict-class
setMethod("[", signature(x = "Rdict"), function(x, i, j, ..., drop = TRUE) {
  for(k in 1:length(x@keys)) {
    if(length(x@keys[[k]]) == length(i) && all(x@keys[[k]] == i))
      return(x@values[[k]])
  }
  stop("key ", i, " was not found")
})

#' @param value The value to assign to key \code{i}.
#' @rdname Rdict-class
setMethod("[<-", signature(x = "Rdict"), function(x, i, j, ..., value) {
  if(is.element(i, x))
    x@values[[i]] <- value
  else {
    x@keys <- c(x@keys, list(i))
    x@values <- c(x@values, list(value))
  }
  return(x)
})

#'
#' The Rdictdefault class.
#'
#' This is a subclass of \linkS4class{Rdict} that contains an additional slot for a default function, which assigns a value to an input key.
#' Only partially implemented, but working well enough for the geometric mean. Will be combined with \linkS4class{Rdict} later.
#'
#' @slot keys A list of keys.
#' @slot values A list of values corresponding to the keys.
#' @slot default A function that takes as input a key and outputs a value to assign to that key.
#' @seealso \linkS4class{Rdict}
#' @name Rdictdefault-class
#' @aliases Rdictdefault
#' @rdname Rdictdefault-class
setClass("Rdictdefault", representation(default = "function"), contains = "Rdict")

#' @param keys A list of keys.
#' @param values A list of values corresponding to the keys.
#' @param default A function that takes as input a key and outputs a value to assign to that key.
#' @rdname Rdictdefault-class
Rdictdefault <- function(keys = list(), values = list(), default) {
  new("Rdictdefault", keys = keys, values = values, default = default)
}

#' @param x A \linkS4class{Rdictdefault} object.
#' @param i A key into the dictionary.
#' @param j,drop,... Unused arguments.
#' @rdname Rdictdefault-class
setMethod("[", signature(x = "Rdictdefault"), function(x, i, j, ..., drop = TRUE) {
  if(length(x@keys) > 0) {
    for(k in 1:length(x@keys)) {
      if(length(x@keys[[k]]) == length(i) && all(x@keys[[k]] == i))
        return(x@values[[k]])
    }
  }
  
  # TODO: Can't update in place. If key doesn't exist, want to create it with default function value.
  stop("Unimplemented: For now, user must manually create key and set its value to default(key)")
  x@keys <- c(x@keys, list(i))
  x@values <- c(x@values, list(x@default(i)))
  return(x@values[[length(x@values)]])
})