Discrete2MixedInt.exprval_in_vec_ineq <- function(expr, vec) {
  expr_dim <- dim(expr)
  if(!(length(expr_dim) == 1 || (length(expr_dim) == 2 && expr_dim[2] == 1)))
    stop("Expression must be 1-dimensional")
  n_entries <- expr_dim[1]

  vec <- sort(vec)
  d <- base::diff(vec)
  repeated_d <- matrix(rep(d, n_entries), nrow = n_entries, byrow = TRUE)
  z <- Variable(nrow(repeated_d), ncol(repeated_d), boolean = TRUE)

  main_con <- (expr == vec[1] + SumEntries(Multiply(repeated_d, z), axis = 1))

  if(size(d) > 1) {
    z_seq <- seq_len(ncol(z) - 1)   # 1, 2, ..., ncol(z) - 1.
    if(length(z_seq) == 0)
      aux_cons <- list()
    else
      aux_cons <- list(z[,z_seq + 1] <= z[,z_seq])
  } else
    aux_cons <- list()
  return(list(main_con, aux_cons))
}

Discrete2MixedInt.exprval_in_vec_eq <- function(expr, vec) {
  expr_dim <- dim(expr)
  if(!(length(expr_dim) == 1 || (length(expr_dim) == 2 && expr_dim[2] == 1)))
    stop("Expression must be 1-dimensional")
  n_entries <- expr_dim[1]

  repeated_vec <- matrix(rep(vec, n_entries), nrow = n_entries, byrow = TRUE)
  z <- Variable(nrow(repeated_vec), ncol(repeated_vec), boolean = TRUE)

  main_con <- (SumEntries(Multiply(repeated_vec, z), axis = 1) == expr)
  aux_cons <- list(SumEntries(z, axis = 1) == 1)
  return(list(main_con, aux_cons))
}

Discrete2MixedInt.get_exprval_in_vec_func <- function(ineq_form) {
  if(ineq_form)
    return(exprval_in_vec_ineq)
  else
    return(exprval_in_vec_eq)
}

Discrete2MixedInt.finite_set_canon <- function(con, .args) {
  vec <- value(con@vec)
  if(size(vec) == 1) {
    # handling for when vec only has a single element
    return(list(con@expre == vec[1], list()))
  }

  flat_expr <- flatten(con@expre)
  exprval_in_vec <- get_exprval_in_vec_func(con@ineq_form)
  res <- exprval_in_vec(flat_expr, vec)
  return(res)
}

#'
#' The Valinvec2MixedInt class.
#'
#' This class represents a discrete to mixed-integer reduction.
#'
#' @rdname Valinvec2MixedInt-class
.Valinvec2MixedInt <- setClass("Valinvec2MixedInt", contains = "Canonicalization")

Valinvec2MixedInt <- function(problem = NULL) { .Valinvec2MixedInt(problem = problem) }

Valinvec2MixedInt.CANON_METHODS <- list("FiniteSet" = Discrete2MixedInt.finite_set_canon)

setMethod("initialize", "Valinvec2MixedInt", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Valinvec2MixedInt.CANON_METHODS)
})

#' @param object A \linkS4class{Valinvec2MixedInt} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Valinvec2MixedInt Checks whether or not the problem involves any finite set constraints.
setMethod("accepts", signature(object = "Valinvec2MixedInt", problem = "Problem"), function(object, problem) {
  return(any(sapply(problem@constraints, function(con) { is(con, "FiniteSet") })))
})
