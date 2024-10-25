## CVXPY SOURCE: cvxpy/reductions/inverse_data.py
#'
#' The InverseData class.
#'
#' This class stores the data useful for solution retrieval.
#'
#' @rdname InverseData-class
.InverseData <- setClass("InverseData", representation(problem = "Problem", id_map = "list", var_offsets = "list", x_length = "numeric", var_dims = "list",
                                                       param_dims = "list", param_to_size = "list", param_id_map = "list",
                                                       id2var = "list", id2cons = "list", cons_id_map = "list", constraints = "ListORNULL"),
                         prototype(id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(),
                                   param_dims = list(), param_to_size = list(), param_id_map = list(), id2var = list(), id2cons = list(),
                                   cons_id_map = list(), constraints = NULL))

InverseData <- function(problem) { .InverseData(problem = problem) }

## Begin R-specific Code Section
## Add InverseData to class union InverseDataORNUL
setIs("InverseData", "InverseDataORNULL")
## End R-specific Code Section

setMethod("initialize", "InverseData", function(.Object, ..., problem, id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(), real2imag = list(), id2cons = list(), cons_id_map = list(), r = NA_real_, minimize = NA, sorted_constraints = list(), is_mip = NA) {
  # Basic variable offset information
  varis <- variables(problem)
  varoffs <- InverseData.get_var_offsets(varis)
  .Object@id_map <- varoffs$id_map
  .Object@var_offsets <- varoffs$var_offsets
  .Object@x_length <- varoffs$vert_offset
  .Object@var_dims <- varoffs$var_dims

  .Object@param_dims <- list()
  # Always start with CONSTANT_ID.
  .Object@param_to_size[[CONSTANT_ID]] <- 1
  .Object@param_id_map <- list()
  offset <- 0
  for(param in parameters(problem)) {
    pid <- as.character(id(param))
    .Object@param_dims[[pid]] <- dim(param)
    .Object@param_to_size[[pid]] <- size(param)
    .Object@param_id_map[[pid]] <- offset
    offset <- offset + size(param)
  }

  # Map of variable id to variable
  .Object@id2var <- stats::setNames(varis, sapply(varis, function(var) { as.character(id(var)) }))

  # Map of constraint id to constraint
  constrs <- problem@constraints
  .Object@id2cons <- stats::setNames(constrs, sapply(constrs, function(cons) { as.character(id(cons)) }))
  .Object@cons_id_map <- list()
  .Object@constraints <- NULL
  return(.Object)
})

InverseData.get_var_offsets <- function(variables) {
  var_dims <- list()
  var_offsets <- list()
  id_map <- list()
  vert_offset <- 0
  for(x in variables) {
    xid <- as.character(id(x))
    var_dims[[xid]] <- dim(x)
    var_offsets[[xid]] <- vert_offset
    id_map[[xid]] <- c(vert_offset, size(x))
    vert_offset <- vert_offset + size(x)
  }
  return(list(id_map = id_map, var_offsets = var_offsets, vert_offset = vert_offset, var_dims = var_dims))
}
