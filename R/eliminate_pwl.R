#'
#' The EliminatePwl class.
#' 
#' This class eliminates piecewise linear atoms.
#' 
#' @rdname EliminatePwl-class
setClass("EliminatePwl", contains = "Canonicalization")

setMethod("accepts", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  atom_types <- sapply(atoms(problem), function(atom) { class(atom) })
  pwl_types <- c("Abs", "MaxElemwise", "SumLargest", "MaxEntries", "Norm1", "NormInf")
  return(any(sapply(atom_types, function(atom) { atom %in% pwl_types })))
})

setMethod("apply", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot canonicalize away piecewise linear atoms.")
  return(apply(Canonicalization(elim_pwl_methods), problem))
})

abs_canon <- function(expr, args) {
  x <- args[[1]]
  t <- Variable(shape(expr))
  constraints <- list(t >= x, t >= -x)
  return(list(t, constraints))
}

max_entries_canon <- function(expr, args) {
  x <- args[[1]]
  shape <- shape(expr)
  axis <- expr@axis
  t <- Variable(shape)
  
  if(is.na(axis))   # shape = c(1,1)
    promoted_t <- promote(t, shape(x))
  else if(axis == 2)   # shape = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = shape(x)[1], ncol = 1) %*% reshape(t, c(1, shape(x)[2])))
  else   # shape = c(m,1)
    promoted_t <- reshape(t, c(shape(x)[1], 1)) %*% Constant(matrix(1, nrow = 1, ncol = shape(x)[2]))
  
  constraints <- list(x <= promoted_t)
  return(list(t, constraints))
}

max_elemwise_canon <- function(expr, args) {
  shape <- shape(expr)
  t <- Variable(shape)
  constraints <- lapply(args, function(elem) { t >= elem })
  return(list(t, constraints))
}

norm1_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  
  # We need an absolute value constraint for the symmetric convex branches (p >= 1)
  constraints <- list()
  
  # TODO: Express this more naturally (recursively) in terms of the other atoms
  abs_expr <- abs(x)
  xconstr <- abs_canon(abs_expr, abs_expr@args)
  abs_x <- xconstr[[1]]
  abs_constraints <- xconstr[[2]]
  constraints <- c(constraints, abs_constraints)
  return(list(SumEntries(abs_x, axis = axis), constraints))
}

norm_inf_canon <- function(expr, args) {
  x <- args[[1]]
  shape <- shape(expr)
  t <- Variable(shape)
  
  promoted_t <- promote(t, shape(x))
  return(list(t, list(x <= promoted_t, x + promoted_t >= 0)))
}

sum_largest_canon <- function(expr, args) {
  x <- args[[1]]
  k <- expr@k
  
  # min sum(t) + kq
  # s.t. x <= t + q, 0 <= t
  t <- Variable(shape(x))
  q <- Variable()
  obj <- sum(t) + k*q
  constraints <- list(x <= t + q, t >= 0)
  return(list(obj, constraints))
}
