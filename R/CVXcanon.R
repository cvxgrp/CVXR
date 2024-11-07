
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
  list(getV = getV, getI = getI, getJ = getJ, getConstVec = getConstVec)
}
