## CVXcanon.ProblemData class shadowing CPP class
        ## problemData = cvxcore.build_matrix(lin_vec,
        ##                                    int(var_length),
        ##                                    id_to_col_C,
        ##                                    param_to_size_C,
        ##                                    s.get_num_threads())

        ##     problemData.param_id = param_id
        ##         problemData.vec_idx = i
        ##         prob_len = problemData.getLen()
        ##         tensor_V[param_id].append(problemData.getV(prob_len))
        ##         tensor_I[param_id].append(problemData.getI(prob_len))
        ##         tensor_J[param_id].append(problemData.getJ(prob_len))

CVXcanon.ProblemData <-R6::R6Class(
  classname = "CVXcanon.ProblemData",
  private = list(
    ptr = NA
  ),
  active = list(
    vec_idx = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_vec_idx', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_vec_idx', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    param_id = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_param_id', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_param_id', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    V = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_V', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_V', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    I = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_I', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_I', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    J = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_J', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_J', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    const_vec = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_const_vec', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_const_vec', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    id_to_col = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_id_to_col', private$ptr, PACKAGE = "CVXR")

      } else {
        .Call('_CVXR_ProblemData__set_id_to_col', private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    const_to_row = function(value) {
      if (missing(value)) {
        .Call('_CVXR_ProblemData__get_const_to_row', private$ptr, PACKAGE = "CVXR")
      } else {
        .Call('_CVXR_ProblemData__set_const_to_row', private$ptr, value, PACKAGE = "CVXR")
      }
    }
  ),
  public = list(
    initialize = function(ptr = NULL) {
      private$ptr <- if (is.null(ptr)) {
        .Call('_CVXR_ProblemData__new', PACKAGE = 'CVXR')
      } else {
        ptr
      }
    }
   ,
    getV = function() {
      .Call('_CVXR_ProblemData__get_V', private$ptr, PACKAGE = "CVXR")
    }
   ,
    getI = function() {
      .Call('_CVXR_ProblemData__get_I', private$ptr, PACKAGE = "CVXR")
    }
   ,
    getJ = function() {
      .Call('_CVXR_ProblemData__get_J', private$ptr, PACKAGE = "CVXR")
    }
   ,
    getConstVec = function() {
      .Call('_CVXR_ProblemData__get_const_vec', private$ptr, PACKAGE = "CVXR")
    }
  ))
