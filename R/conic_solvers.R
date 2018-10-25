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
