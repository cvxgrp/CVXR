## CVXPY SOURCE: cvxpy/utilities/debug_tools.py
#######################
#                     #
# Debugging utilities #
#                     #
#######################
node_count <- function(expr) {
  # Return node count for the expression/constraint.
  if("args" %in% slotNames(expr) && length(expr@args) >= 1) {
    return(1 + sum(sapply(expr@args, node_count)))
  } else
    return(1)
}

build_non_disciplined_error_msg <- function(problem, discipline_type) {
  prop_name <- NA_character_
  prefix_conv <- ""
  if(discipline_type == "DCP")
    prop_name <- "is_dcp"
  else if(discipline_type == "DGP") {
    prop_name <- "is_dgp"
    prefix_conv <- "log_log_"
  } else
    stop("Unknown discipline type ", discipline_type)

  find_non_prop_leaves <- function(expr, res = NULL) {
    if(is.null(res))
      res <- c()
    # if((is.null(expr@args) || length(expr@args) == 0) && slot(expr, prop_name)())
    if((is.null(expr@args) || length(expr@args) == 0) && do.call(prop_name, expr))
      return(res)

    # if(!slot(expr, prop_name)() && all(sapply(expr@args, function(child) { slot(child, prop_name)() } ))) {
    if(!do.call(prop_name, expr) && all(sapply(expr@args, function(child) { do.call(prop_name, child) } ))) {
      str_expr <- as.character(expr)
      if(discipline_type == "DGP" && is(expr, "Variable"))
        str_expr <- paste(str_expr, "<-- needs to be declared positive")
      res <- c(res, str_expr)
    }

    for(child in expr@args)
      res <- find_non_prop_leaves(child, res)
    return(res)
  }

  # if(!slot(problem@objective, prop_name)()) {
  if(!do.call(prop_name, problem@objective)) {
    non_disciplined_leaves <- find_non_prop_leaves(problem@objective@expr)
    if(length(non_disciplined_leaves) > 0)
      msg <- paste("The objective is not ", discipline_type, ". Its following subexpressions are not:", sep = "")
    else {
      convex_str <- paste(prefix_conv, "convex", sep = "")
      concave_str <- paste(prefix_conv, "concave", sep = "")
      # fun_attr_check <- slot(problem@objective@args[[1]], paste("is_", convex_str, sep = ""))()
      fun_attr_check <- do.call(paste("is_", convex_str, sep = ""), problem@objective@args[[1]])
      msg <- paste("The objective is not ", discipline_type, ", even though each sub-expression is.\n",
                   "You are trying to ", problem@objective@NAME, " a function that is ",
                   ifelse(fun_attr_check, convex_str, concave_str), sep = "")
    }

    for(expr in non_disciplined_leaves)
      msg <- paste(msg, as.character(expr), sep = "\n")
    return(msg)
  }

  disciplined_mask <- sapply(problem@constraints, is_dcp)
  not_disciplined_constraints <- problem@constraints[!disciplined_mask]

  msg <- paste("The following constraints are not", discipline_type)
  for(expr in not_disciplined_constraints) {
    msg <- paste(msg, "\n", as.character(expr), ", because the following subexpressions are not:", sep = "")
    non_disciplined_leaves <- find_non_prop_leaves(expr)
    for(subexpr in non_disciplined_leaves)
      msg <- paste(msg, "\n|-- ", as.character(subexpr), sep = "")
  }
  return(msg)
}
