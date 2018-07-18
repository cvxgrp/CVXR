# Factory function for infeasible or unbounded solutions.
failure_solution <- function(status) {
  if(status == INFEASIBLE)
    opt_val = Inf
  else if(status == UNBOUNDED)
    opt_val = -Inf
  else
    opt_val = NA_real_
  return(Solution(status, opt_val, list(), list(), list()))
}

#'
#' The Solution class.
#'
#' This class represents a solution to an optimization problem.
#' 
#' @rdname Solution-class
.Solution <- setClass("Solution", representation(status = "character", opt_val = "numeric", primal_vars = "list", dual_vars = "list", attr = "list"),
                     prototype(primal_vars = list(), dual_vars = list(), attr = list()))

Solution <- function(status, opt_val, primal_vars, dual_vars, attr) {
  .Solution(status = status, opt_val = opt_val, primal_vars = primal_vars, dual_vars = dual_vars, attr = attr)
}

setMethod("show", "Solution", function(object) {
  cat("Solution(", object@status, ", (", 
                  paste(object@primal_vars, collapse = ", "), "), (", 
                  paste(object@dual_vars, collapse = ", "), "), (", 
                  paste(object@attr, collapse = ", "), "))", sep = "")
})

setMethod("as.character", "Solution", function(x) {
  paste("Solution(", x@status, ", (", 
                    paste(x@primal_vars, collapse = ", "), "), (", 
                    paste(x@dual_vars, collapse = ", "), "), (", 
                    paste(x@attr, collapse = ", "), "))", sep = "")
})

#'
#' The InverseData class.
#'
#' This class represents the data encoding an optimization problem.
#' 
#' @rdname InverseData-class
.InverseData <- setClass("InverseData", representation(problem = "Problem", id_map = "list", var_offsets = "list", x_length = "numeric", var_shapes = "list",
                                                       id2var = "list", real2imag = "list", id2cons = "list", cons_id_map = "list"),
                                        prototype(id_map = list(), var_offsets = list(), x_length = NA_real_, var_shapes = list(), id2var = list(),
                                                  real2imag = list(), id2cons = list(), cons_id_map = list()))

InverseData <- function(problem) { .InverseData(problem = problem) }

setMethod("initialize", "InverseData", function(.Object, ..., problem, id_map = list(), var_offsets = list(), x_length = NA_real_, var_shapes = list(), id2var = list(), real2imag = list(), id2cons = list(), cons_id_map = list()) {
  # Basic variable offset information
  varis <- variables(problem)
  varoffs <- get_var_offsets(.Object, varis)
  .Object@id_map <- varoffs$id_map
  .Object@var_offsets <- varoffs$var_offsets
  .Object@x_length <- varoffs$x_length
  .Object@var_shapes <- varoffs$var_shapes
  
  # Map of variable id to variable
  .Object@id2var <- setNames(varis, sapply(varis, function(var) { as.character(id(var)) }))
  
  # Map of real to imaginary parts of complex variables
  var_comp <- lapply(varis, function(var) { if(is_complex(var)) var })
  .Object@real2imag <- setNames(var_comp, sapply(var_comp, function(var) { as.character(id(var)) }))
  constrs <- constraints(problem)
  constr_comp <- lapply(constrs, function(cons) { if(is_complex(cons)) cons })
  constr_dict <- setNames(constr_comp, sapply(constr_comp, function(cons) { as.character(id(cons)) }))
  .Object@real2imag <- update(.Object@real2imag, constr_dict)
  
  # Map of constraint id to constraint
  .Object@id2cons <- setNames(constrs, sapply(constrs, function(cons) { id(cons) }))
  .Object@cons_id_map <- list()
})

setMethod("get_var_offsets", signature(object = "InverseData", variables = "list"), function(object, variables) {
  var_shapes <- list()
  var_offsets <- list()
  id_map <- list()
  vert_offset <- 0
  for(x in variables) {
    var_shapes[[as.character(id(x))]] <- shape(x)   # TODO: Redefine Variable class to include shape parameter
    var_offsets[[as.character(id(x))]] <- vert_offset
    id_map[[as.character(id(x))]] <- list(vert_offset, size(x))
    vert_offset <- vert_offset + size(x)
  }
  return(list(id_map = id_map, var_offsets = var_offsets, x_length = vert_offset, var_shapes = var_shapes))
})

#'
#' The Reduction class.
#'
#' This virtual class represents a reduction, an actor that transforms a problem
#' into an equivalent problem. By equivalent, we mean that there exists a mapping
#' between solutions of either problem: if we reduce a problem \eqn{A} to another
#' problem \eqn{B} and then proceed to find a solution to \eqn{B}, we can convert
#' it to a solution of \eqn{A} with at most a moderate amount of effort.
#' 
#' Every reduction supports three methods: accepts, apply, and invert. The accepts
#' method of a particular reduction codifies the types of problems that it is applicable
#' to, the apply method takes a problem and reduces it to a (new) equivalent form,
#' and the invert method maps solutions from reduced-to problems to their problems
#' of provenance.
#' 
#' @rdname Reduction-class
setClass("Reduction", contains = "VIRTUAL")

#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Reduction States whether the reduction accepts a problem.
setMethod("accepts", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @describeIn Reduction Applies the reduction to a problem and returns an equivalent problem.
setMethod("apply", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The data encoding the original problem.
#' @describeIn Reduction Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Reduction", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { stop("Unimplemented") })
