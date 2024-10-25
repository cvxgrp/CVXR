## CVXPY SOURCE: cvxpy/utilities/versioning.py
#'
#' The Version class.
#'
#' The class represents the version of CVXR.
#' 
#' It parses strings of the form 'x.y.z+[STUFF]' and tuples of the form (x, y, z).
#' In the former case, '+[STUFF]' refers to a local version identifier in the
#' sense of https://www.python.org/dev/peps/pep-0440/#local-version-identifiers.
#' We don't actually store the local version identifier for comparisons.
#'
#' @name Version-class
#' @aliases Version
#' @rdname Version-class
.Version <- setClass("Version", representation(v = "ANY", major = "integer", minor = "integer", micro = "integer"),
                                prototype(major = NA_integer_, minor = NA_integer_, micro = NA_integer_))

Version <- function(v) { new("Version", v = v) }

setMethod("initialize", "Version", function(.Object, ..., v, major, minor, micro) {
  if(is.character(v)) {
    v <- strsplit(v, "rc")[[1]][1]
    v <- strsplit(v, ".")[[1]]
    if(length(v) < 3)
      stop("v must be splittable into at least 3 places")   # Anything after the third place doesn't matter.
    v[3] <- strsplit(v[3], "+")[[1]][1]   # Anything after the + doesn't matter.
  }
  
  .Object@major <- as.integer(v[1])
  .Object@minor <- as.integer(v[2])
  .Object@micro <- as.integer(v[3])
  .Object@v <- c(.Object@major, .Object@minor, .Object@micro)
  return(.Object)
})

# Comparison operators
#' @param e1,e2 The \linkS4class{Version} objects or numeric constants to compare.
#' @rdname Version-class
setMethod("==", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v == e2@v })

#' @rdname Version-class
setMethod("!=", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v != e2@v })

#' @rdname Version-class
setMethod("<=", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v <= e2@v })

#' @rdname Version-class
setMethod("<", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v < e2@v })

#' @rdname Version-class
setMethod(">=", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v >= e2@v })

#' @rdname Version-class
setMethod(">", signature(e1 = "Version", e2 = "Version"), function(e1, e2) { e1@v > e2@v })

#' @describeIn Version The string representation of the version.
setMethod("as.character", "Version", function(x) {
  paste(x@v, collapse = ".")
})
