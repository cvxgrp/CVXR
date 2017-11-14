## CVXcanon.ProblemData class shadowing CPP class
CVXcanon.ProblemData <- R6::R6Class("CVXcanon.ProblemData",
                                    private = list(
                                        ptr = NA
                                    ),
                                    active = list(
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
