## An interface to build_matrix in C
## Return a ProblemData object
CVXcanon.build_matrix <- function(constraints, id_to_col, constr_offsets) {
  ## Returns a ProblemData object
  ptr <- NA
  if (missing(constr_offsets)) {
    ptr <- .Call('_CVXR_build_matrix_0',
                 constraints$getXPtr(),
                 id_to_col,
                 PACKAGE = 'CVXR')
  } else {
    ptr <- .Call('_CVXR_build_matrix_1',
                 constraints$getXPtr(),
                 id_to_col,
                 constr_offsets,
                 PACKAGE = 'CVXR')
  }

  getV <- function() {
    .Call('_CVXR_ProblemData__get_V', ptr, PACKAGE = "CVXR")
  }

  getI <- function() {
    .Call('_CVXR_ProblemData__get_I', ptr, PACKAGE = "CVXR")
  }

  getJ <- function() {
    .Call('_CVXR_ProblemData__get_J', ptr, PACKAGE = "CVXR")
  }

  getConstVec <- function() {
    .Call('_CVXR_ProblemData__get_const_vec', ptr, PACKAGE = "CVXR")
  }
  
  getLen <- function() {
    .Call('_CVXR_ProblemData__get_len', ptr, PACKAGE = "CVXR")
  }
    
  set_vec_idx <- function(int_value) {
    .Call('_CVXR_ProblemData__set_vec_idx', ptr, int_value, PACKAGE = "CVXR")    
  }

  get_vec_idx <- function() {
    .Call('_CVXR_ProblemData__get_vec_idx', ptr, PACKAGE = "CVXR")    
  }

  set_param_id <- function(int_value) {
    .Call('_CVXR_ProblemData__set_param_id', ptr, int_value, PACKAGE = "CVXR")    
  }

  get_param_id <- function() {
    .Call('_CVXR_ProblemData__get_param_id', ptr, PACKAGE = "CVXR")    
  }

  ## Return a ProblemData object
  list(getV = getV, getI = getI, getJ = getJ, getConstVec = getConstVec, getLen = getLen,
       set_vec_idx = set_vec_idx, get_vec_idx = get_vec_idx,
       set_param_id = set_param_id, get_param_id = get_param_id)
}
