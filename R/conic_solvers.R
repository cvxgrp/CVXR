is_stuffed_cone_constraint <- function(constraint) {
  # Conic solvers require constraints to be stuffed in the following way.
  if(length(variables(constraint)) != 1)
    return(FALSE)
  for(arg in constraints@args) {
    if(class(arg) == "Reshape")
      arg <- arg@args[[1]]
    if(class(arg) == "AddExpression") {
      if(!class(arg@args[[1]]) %in% c("MulExpression", "Multiply"))
        return(FALSE)
      if(class(arg@args[[1]]@args[[1]]) != "Constant")
        return(FALSE)
      if(class(arg@args[[2]]) != "Constant")
        return(FALSE)
    } else if(class(arg) %in% c("MulExpression", "Multiply")) {
      if(class(arg@args[[1]]) != "Constant")
        return(FALSE)
    } else
      return(FALSE)
  }
  return(TRUE)
}

is_stuffed_cone_objective <- function(objective) {
  # Conic solvers require objectives to be stuffed in the following way.
  expr <- objective@expr
  return(is_affine(expr) && length(variables(expr)) == 1 && class(expr) == "AddExpression" && length(expr@args) == 2
                         && class(expr@args[[1]]) %in% c("MulExpression", "Multiply") && class(expr@args[[2]]) == "Constant")
}


#' Summary of cone dimensions present in constraints.
#'
#'    Constraints must be formatted as dictionary that maps from
#'    constraint type to a list of constraints of that type.
#'
#'    Attributes
#'    ----------
#'    zero : int
#'        The dimension of the zero cone.
#'    nonpos : int
#'        The dimension of the non-positive cone.
#'    exp : int
#'        The dimension of the exponential cone.
#'    soc : list of int
#'        A list of the second-order cone dimensions.
#'    psd : list of int
#'        A list of the positive semidefinite cone dimensions, where the
#'        dimension of the PSD cone of k by k matrices is k.
.ConeDims <- setClass("ConeDims", representation(constr_map = "list", zero = "numeric", nonpos = "numeric", exp = "numeric", soc = "numeric", psd = "numeric"),
                                  prototype(zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = NA_real_, psd = NA_real_))

ConeDims <- function(constr_map) { .ConeDims(constr_map = constr_map) }

setMethod("initialize", "ConeDims", function(.Object, constr_map, zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = NA_real_, psd = NA_real_) {
  .Object@zero <- sum(sapply(constr_map$Zero, function(c) { size(c) }))
  .Object@nonpos <- sum(sapply(constr_map$NonPos, function(c) { size(c) }))
  .Object@exp <- sum(sapply(constr_map$ExpCone, function(c) { num_cones(c) }))
  .Object@soc <- sapply(constr_map$SOC, function(c) { cone_sizes(c) })
  .Object@psd <- sapply(constr_map$PSD, function(c) { dim(c)[1] })
  return(.Object)
})

# Conic solver class with reduction semantics.
ConicSolver <- setClass("ConicSolver", contains = "Solver")

# The key that maps to ConeDims in the data returned by apply().
setMethod("dims", "ConicSolver", function(object) { "dims" })

# Every conic solver must support Zero and NonPos constraints.
setMethod("supported_constraints", "ConicSolver", function(object) { c("Zero", "NonPos") })

# Some solvers cannot solve problems that do not have constraints.
# For such solvers, requires_constr should return TRUE.
setMethod("requires_constr", "ConicSolver", function(object) { FALSE })

setMethod("accepts", signature(object = "ConicSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && (mip_capable(object) || !is_mixed_integer(problem)) && is_stuffed_cone_objective(problem@objective)
    && !convex_attributes(variables(problem)) && (length(problem@constraints) > 0 || !requires_constr(object))
    && all(sapply(problem@constraints, function(c) { class(c) %in% supported_constraints(object) }))
    && all(sapply(problem@constraints, function(c) { is_stuffed_cone_constraints(c) })))
})

ConicSolver.get_coeff_offset(expr) {
  if(class(expr) == "Reshape")
    expr <- expr@args[[1]]
  if(length(expr@args[[1]]@args) == 0) {
    offset <- 0
    coeff <- as.numeric(value(expr@args[[1]]))
  } else {
    offset <- as.numeric(value(expr@args[[2]]))
    coeff <- as.numeric(value(expr@args[[1]]@args[[1]]))
  }
  # TODO: Finish this implementation and port rest of ConicSolver.
}
