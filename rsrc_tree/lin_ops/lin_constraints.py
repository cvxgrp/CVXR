## CVXPY SOURCE: cvxpy/lin_ops/lin_constraints.py

LinConstr <- function(expr, constr_id, class = "LinConstr") {
  if (is.null(dim)) dim <- c(1L, 1L)   # TODO: Get rid of this with proper multi-dimensional handling.
  ##if(!is.character(constr_id)) stop("constr_id must be a character string")
  if (!is.integer(constr_id)) stop("constr_id must be an integer")
  if (!is.numeric(dim)) stop("dim must be a numeric vector")
  list(expr = expr, constr_id = constr_id, dim = dim, class = class)
}

LinEqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinEqConstr") }
LinLeqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinLeqConstr") }

