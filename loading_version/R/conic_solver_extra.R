#'
#' Is the constraint a stuffed cone constraint?
#'
#' @param constraint A \linkS4class{Constraint} object.
#' @return Is the constraint a stuffed-cone constraint?
is_stuffed_cone_constraint <- function(constraint) {
  # Conic solvers require constraints to be stuffed in the following way.
  if(length(variables(constraint)) != 1)
    return(FALSE)
  for(arg in constraint@args) {
    if(inherits(arg, "Reshape"))
      arg <- arg@args[[1]]
    if(inherits(arg, "AddExpression")) {
      if(!inherits(arg@args[[1]], c("MulExpression", "Multiply")))
        return(FALSE)
      if(!inherits(arg@args[[1]]@args[[1]], "Constant"))
        return(FALSE)
      if(!inherits(arg@args[[2]], "Constant"))
        return(FALSE)
    } else if(inherits(arg, c("MulExpression", "Multiply"))) {
      if(!inherits(arg@args[[1]], "Constant"))
        return(FALSE)
    } else
      return(FALSE)
  }
  return(TRUE)
}

#'
#' Is the objective a stuffed cone objective?
#'
#' @param objective An \linkS4class{Objective} object.
#' @return Is the objective a stuffed-cone objective?
is_stuffed_cone_objective <- function(objective) {
  # Conic solvers require objectives to be stuffed in the following way.
  expr <- expr(objective)
  return(is_affine(expr) && length(variables(expr)) == 1 && inherits(expr, "AddExpression") && length(expr@args) == 2
                         && inherits(expr@args[[1]], c("MulExpression", "Multiply")) && inherits(expr@args[[2]], "Constant"))
}

