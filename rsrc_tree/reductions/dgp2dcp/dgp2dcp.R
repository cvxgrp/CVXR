## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/dgp2dcp.py
#'
#' Reduce DGP problems to DCP problems.
#'
#' This reduction takes as input a DGP problem and returns an equivalent DCP
#' problem. Because every (generalized) geometric program is a DGP problem,
#' this reduction can be used to convert geometric programs into convex form.
#' @rdname Dgp2Dcp-class
.Dgp2Dcp <- setClass("Dgp2Dcp", contains = "Canonicalization")
Dgp2Dcp <- function(problem = NULL) { .Dgp2Dcp(problem = problem) }

setMethod("initialize", "Dgp2Dcp", function(.Object, ..., problem = NULL) {
  # Canonicalization of DGP is stateful; canon_methods is created in 'perform'.
  callNextMethod(.Object, ..., canon_methods = NULL, problem = problem)
})

#' @param object A \linkS4class{Dgp2Dcp} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Dgp2Dcp Is the problem DGP?
setMethod("accepts", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  is_dgp(problem) &&
  all(sapply(parameters(problem), function(p) {
      p_val <- value(p)
      return(!is.null(p_val) && !is.na(p_val))
    }))
})

#' @describeIn Dgp2Dcp Converts the DGP problem to a DCP problem.
setMethod("perform", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("The supplied problem is not DGP")

  object@canon_methods <- DgpCanonMethods()
  tmp <- callNextMethod(object, problem)
  object <- tmp[[1]]
  equiv_problem <- tmp[[2]]
  inverse_data <- tmp[[3]]
  inverse_data@problem <- problem
  return(list(object, equiv_problem, inverse_data))
})

#' @param expr An \linkS4class{Expression} object corresponding to the DGP problem.
#' @param args A list of values corresponding to the DGP expression
#' @describeIn Dgp2Dcp Canonicalizes each atom within an Dgp2Dcp expression.
setMethod("canonicalize_expr", "Dgp2Dcp", function(object, expr, args) {
  if(inherits(expr, names(object@canon_methods)))
    return(object@canon_methods[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn Dgp2Dcp Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "Dgp2Dcp", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  solution <- callNextMethod(object, solution, inverse_data)
  if(solution@status == SOLVER_ERROR)
    return(solution)
  for(vid in names(solution@primal_vars))
    solution@primal_vars[[vid]] <- exp(solution@primal_vars[[vid]])
  # f(x) = e^{F(u)}
  solution@opt_val <- exp(solution@opt_val)
  return(solution)
})
