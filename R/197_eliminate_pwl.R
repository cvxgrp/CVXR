## CVXPY SOURCE: cvxpy/reductions/eliminate_pwl.py
#'
#' The EliminatePwl class.
#'
#' This class eliminates piecewise linear atoms.
#'
#' @rdname EliminatePwl-class
.EliminatePwl <- setClass("EliminatePwl", contains = "Canonicalization")
EliminatePwl <- function(problem = NULL) { .EliminatePwl(problem = problem) }

setMethod("initialize", "EliminatePwl", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = EliminatePwl.CANON_METHODS)
})

#' @param object An \linkS4class{EliminatePwl} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn EliminatePwl Does this problem contain piecewise linear atoms?
setMethod("accepts", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  atom_types <- sapply(atoms(problem), function(atom) { class(atom) })
  pwl_types <- c("Abs", "MaxElemwise", "SumLargest", "MaxEntries", "Norm1", "NormInf")
  return(any(sapply(atom_types, function(atom) { atom %in% pwl_types })))
})

setMethod("perform", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot canonicalize away piecewise linear atoms.")
  callNextMethod(object, problem)
})

