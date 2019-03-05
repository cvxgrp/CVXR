#'
#' The EliminatePwl class.
#' 
#' This class eliminates piecewise linear atoms.
#' 
#' @rdname EliminatePwl-class
.EliminatePwl <- setClass("EliminatePwl", representation(problem = "Problem"), prototype(problem = NULL), contains = "Canonicalization")

EliminatePwl <- function(problem = NULL) { .EliminatePwl(problem = problem) }

setMethod("initialize", "EliminatePwl", function(.Object, ..., problem = NULL) {
  .Object@problem <- problem
  callNextMethod(.Object, ..., problem = problem, canon_methods = EPWL_CANON_METHODS)
})

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

# Atom canonicalizers.
EliminatePwl.abs_canon <- function(expr, args) {
  x <- args[[1]]
  t <- Variable(dim(expr))
  constraints <- list(t >= x, t >= -x)
  return(list(t, constraints))
}

EliminatePwl.max_entries_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  axis <- expr@axis
  t <- Variable(expr_dim)
  
  if(is.na(axis))   # dim(expr) = c(1,1)
    promoted_t <- promote(t, dim(x))
  else if(axis == 2)   # dim(expr) = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = 1) %*% reshape(t, c(1, ncol(x))))
  else   # shape = c(m,1)
    promoted_t <- reshape(t, c(nrow(x), 1)) %*% Constant(matrix(1, nrow = 1, ncol = ncol(x)))
  
  constraints <- list(x <= promoted_t)
  return(list(t, constraints))
}

EliminatePwl.max_elemwise_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)
  constraints <- lapply(args, function(elem) { t >= elem })
  return(list(t, constraints))
}

EliminatePwl.min_entries_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Length of args must be one")
  tmp <- MaxEntries(-args[[1]])
  canon <- EliminatePwl.max_entries_canon(tmp, tmp@args)
  return(list(-canon[[1]], canon[[2]]))
}

EliminatePwl.min_elemwise_canon <- function(expr, args) {
  tmp <- do.call(MaxElemwise, lapply(args, function(arg) { -arg }))
  canon <- EliminatePwl.max_elemwise_canon(tmp, tmp@args)
  return(list(-canon[[1]], canon[[2]]))
}

EliminatePwl.norm1_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  
  # We need an absolute value constraint for the symmetric convex branches (p >= 1)
  constraints <- list()
  
  # TODO: Express this more naturally (recursively) in terms of the other atoms
  abs_expr <- abs(x)
  xconstr <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_x <- xconstr[[1]]
  abs_constraints <- xconstr[[2]]
  constraints <- c(constraints, abs_constraints)
  return(list(SumEntries(abs_x, axis = axis), constraints))
}

EliminatePwl.norm_inf_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)
  
  promoted_t <- promote(t, dim(x))
  return(list(t, list(x <= promoted_t, x + promoted_t >= 0)))
}

EliminatePwl.sum_largest_canon <- function(expr, args) {
  x <- args[[1]]
  k <- expr@k
  
  # min sum(t) + kq
  # s.t. x <= t + q, 0 <= t
  t <- Variable(dim(x))
  q <- Variable()
  obj <- sum(t) + k*q
  constraints <- list(x <= t + q, t >= 0)
  return(list(obj, constraints))
}

EliminatePwl.CANON_METHODS <- list(Abs = EliminatePwl.abs_canon,
                                   MaxElemwise = EliminatePwl.max_elemwise_canon,
                                   MaxEntries = EliminatePwl.max_entries_canon,
                                   MinElemwise = EliminatePwl.min_elemwise_canon,
                                   MinEntries = EliminatePwl.min_entries_canon,
                                   Norm1 = EliminatePwl.norm1_canon,
                                   NormInf = EliminatePwl.norm_inf_canon,
                                   SumLargest = EliminatePwl.sum_largest_canon)
