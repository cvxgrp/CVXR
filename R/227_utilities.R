## CVXPY SOURCE: cvxpy/reductions/utilities.py

#==============================#
# Reduction utility functions
#==============================#
lower_ineq_to_nonpos <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonPosConstraint(lhs - rhs, constr_id = inequality@constr_id)
}

lower_ineq_to_nonneg <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonNegConstraint(rhs - lhs, constr_id = inequality@constr_id)
}

lower_equality <- function(equality) {
  lhs <- equality@args[[1]]
  rhs <- equality@args[[2]]
  ZeroConstraint(lhs - rhs, constr_id = equality@constr_id)
}

nonpos2nonneg <- function(nonpos) {
  NonNegConstraint(-nonpos@expr, constr_id = nonpos@constr_id)
}

special_index_canon <- function(expr, args) {
  select_mat <- expr@.select_mat
  final_dim <- dim(expr@.select_mat)
  if(is.null(final_dim))
    final_dim <- c(length(select_mat), 1)
  # select_vec <- matrix(select_mat, nrow = prod(final_dim))
  select_vec <- as.vector(select_mat)

  # Select the chosen entries from expr.
  arg <- args[[1]]
  arg_size <- size(arg)
  # identity <- diag(size(arg))
  identity <- sparseMatrix(i = seq(arg_size), j = seq(arg_size), x = rep(1, arg_size))

  v <- vec(arg)
  idmat <- matrix(identity[select_vec,], ncol = arg_size)
  if(is_scalar(v) || is_scalar(as.Constant(idmat)))
    lowered <- Reshape(idmat * v, final_dim)
  else
    lowered <- Reshape(idmat %*% v, final_dim)
  list(lowered, list())
}

#'
#' Are the arguments affine?
#'
#' @param constraints A list of \linkS4class{Constraint} object.
#' @return All the affine arguments in given constraints.
are_args_affine <- function(constraints) {
  all(sapply(constraints, function(constr) { all(sapply(constr@args, is_affine)) }))
}


## #' REPLACED BELOW by modified version
## #' Organize the constraints into a list keyed by constraint names.
## #'
## #'@param constraints A list of \linkS4class{Constraint} objects
## #'@return A list keyed by constraint types where list[[cone_type]]
## #'   maps to a list of exactly those constraints that are of type
## #'   cone_type.
## group_constraints <- function(constraints) {
##   constr_map <- list()
##   for(c in constraints) {
##     if(class(c) %in% names(constr_map))
##       constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
##     else
##       constr_map[[class(c)]] <- list(c)
##   }
##   return(constr_map)
## }


#'
#' Organize the constraints into a list keyed by constraint names.
#'
#'@param constraints A list of \linkS4class{Constraint} objects
#'@return A list keyed by constraint types where list[[cone_type]] maps to a list of exactly those constraints that are of type cone_type.
group_constraints <- function(constraints) {
  ## The constr_types list below should match the map named used in ConeDims-class (file dcp2cone.R)
  constr_map <- list()
  constr_names <- character(0)
  for (constr in constraints) {
    cl <- class(constr)
    index <- match(cl, constr_names)
    if (is.na(index)) { ## class not yet appeared in list
      constr_map[[cl]] <- list(constr) ## add new named item to list
      constr_names <- c(constr_names, cl)  ## add name to list of names
    } else {
      constr_map[[cl]] <- c(constr_map[[cl]], list(constr))
    }
  }
  constr_map
}

#'
#' The ReducedMat class.
#'
#' This is a utility class for condensing the mapping from parameters to problem data.
#'
#' For maximum efficiency of representation and application, the mapping from
#' parameters to problem data must be condensed. It begins as a CSC sparse matrix
#' matrix_data, such that multiplying by a parameter vector gives the problem data.
#' The row index array and column pointer array are saved as problem_data_index,
#' and a CSR matrix reduced_mat that when multiplied by a parameter vector gives
#' the values array. The ReducedMat class caches the condensed representation
#' and provides a method for multiplying by a parameter vector.
#'
#' This class consolidates code from ParamConeProg and ParamQuadProg.
#'
#' @rdname ReducedMat-class
ReducedMat <- setClass("ReducedMat", representation(matrix_data = "numeric", var_len = "integer", quad_form = "logical", reduced_mat = "numeric", problem_data_index = "ListORNULL", mapping_nonzero = "numeric"),
                        prototype(quad_form = FALSE, reduced_mat = NA_real_, problem_data_index = NULL, mapping_nonzero = NA_integer_))

#' @param keep_zeros (Optional) If TRUE, store explicit zeros in A where parameters are affected.
#' @describeIn ReducedMat Cache computed attributes if not present.
setMethod("cache", "ReducedMat", function(object, keep_zeros = FALSE) {
  # Short circuit null case.
  if(is.na(object@matrix_data))
    return(object)

  if(is.na(object@reduced_mat)) {
    # Form a reduced representation of the mapping, for faster application of parameters.
    if(!is.null(dim(object@matrix_data)) && prod(dim(object@matrix_data)) != 0) {
      tmp <- canonInterface.reduce_problem_data_tensor(object@matrix_data, object@var_len, object@quad_form)
      reduced_mat <- tmp[[1]]
      indices <- tmp[[2]]
      indptr <- tmp[[3]]
      shape <- tmp[[4]]
      object@reduced_mat <- reduced_mat
      object@problem_data_index <- list(indices, indptr, shape)
    } else {
      object@reduced_mat <- object@matrix_data
      object@problem_data_index <- NULL
    }
  }

  if(keep_zeros && is.na(object@mapping_nonzero)) {
    object@mapping_nonzero <- canonInterface.A_mapping_nonzero_rows(object@matrix_data, object@var_len)
  }
  return(object)
})

#' @param param_vec Flattened parameter vector
#' @param with_offset (Optional) A logical value indicating whether to return offset. Defaults to TRUE.
#' @describeIn ReducedMat Wraps get_matrix_from_tensor in canonInterface
setMethod("get_matrix_from_tensor", "ReducedMat", function(object, param_vec, with_offset = TRUE) {
  canonInterface.get_matrix_from_tensor(object@reduced_mat, param_vec, object@var_len, nonzero_rows = object@mapping_nonzero, with_offset = with_offset, problem_data_index = object@problem_data_index)
})

## Begin R-specific Code
setIs("ReducedMat", "ReducedMatORNULL")
## End R-specific Code
