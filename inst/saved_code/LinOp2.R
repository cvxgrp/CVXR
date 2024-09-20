## Another implementation
#' An R LinOp class as distinguished from a C++ class of the same name. Always has a unique id
#' @param type the type of LinOp, one of the types above
#' @param dim the shape of the LinOp, a tuple, so for us a vector of integers
#' @param args the arguments of the LinOp
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @importFrom uuid UUIDgenerate
#' @return an object of class "LinOp"
LinOp <- function(type, dim, args, data = NULL, id = uuid::UUIDgenerate()) {
  get_type <- function() type
  set_type <- function(what) type <<- what
  get_dim <- function() dim
  set_dim <- function(what) dim <<- what
  get_args <- function() args
  set_args <- function(what) args <<- what
  result <- list(get_type = get_type, set_type = set_type,
                 get_dim = get_dim, set_dim = set_dim,
                 get_args = get_args, set_args = set_args,
                 get_id = function() id)
  result$self <- result
  class(result) <- "LinOp"
  result
}

#' @method print LinOp
print.LinOp <- function(x, ...) {
  sprintf("LinOp(%s, dim = [%s])",
          x$self$get_type(),
          paste0(x$self$get_dim(), collapse = ", ")
          )
}

#' Make a (R) Linear Constraint
#' @param expr the expression
#' @param constr_id the constaint id
#' @param dim the shape

#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @return an object of class "LinOp"
make_lin_constraint <- function(expr, constr_id, dim, class = "LinConstr") {
  result <- list(get_expr = function() expr
                 get_constr_id = function() constr_id,
                 get_dim = function() dim)
  result$self <- result
  class(result) <- class
  result
}

LinConstr <- make_lin_constraint
LinEqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinEqConstr") }
LinLeqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinLeqConstr") }

print_lin_constraint <- function(x, ...) {
  sprintf("%s(%s, dim = [%s])",
          class(x),
          x$self$get_expr(),
          paste0(x$self$get_dim(), collapse = ", ")
          )
}

#' @method print LinConstr
print.LinConstr <- print_lin_constraint

#' @method print LinEqConstr
print.LinEqConstr <- print_lin_constraint

#' @method print LinLeqConstr
print.LinLeqConstr <- print_lin_constraint






