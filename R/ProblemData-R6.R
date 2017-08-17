## CVXcanon.ProblemData class shadowing CPP class
CVXcanon.ProblemData <- R6::R6Class("CVXcanon.ProblemData",
                                    private = list(
                                        ptr = NA
                                    ),
                                    active = list(
                                        V = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_V', private$ptr, PACKAGE = "cvxr")
                                            } else {
                                                .Call('_cvxr_ProblemData__set_V', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                       ,
                                        I = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_I', private$ptr, PACKAGE = "cvxr")
                                            } else {
                                                .Call('_cvxr_ProblemData__set_I', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                       ,
                                        J = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_J', private$ptr, PACKAGE = "cvxr")
                                            } else {
                                                .Call('_cvxr_ProblemData__set_J', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                       ,
                                        const_vec = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_const_vec', private$ptr, PACKAGE = "cvxr")
                                            } else {
                                                .Call('_cvxr_ProblemData__set_const_vec', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                       ,
                                        id_to_col = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_id_to_col', private$ptr, PACKAGE = "cvxr")

                                            } else {
                                                .Call('_cvxr_ProblemData__set_id_to_col', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                       ,
                                        const_to_row = function(value) {
                                            if (missing(value)) {
                                                .Call('_cvxr_ProblemData__get_const_to_row', private$ptr, PACKAGE = "cvxr")
                                            } else {
                                                .Call('_cvxr_ProblemData__set_const_to_row', private$ptr, value, PACKAGE = "cvxr")
                                            }
                                        }
                                    ),
                                    public = list(
                                        initialize = function(ptr = NULL) {
                                            private$ptr <- if (is.null(ptr)) {
                                                               .Call('_cvxr_ProblemData__new', PACKAGE = 'cvxr')
                                                           } else {
                                                               ptr
                                                           }
                                        }
                                       ,
                                        getV = function() {
                                                .Call('_cvxr_ProblemData__get_V', private$ptr, PACKAGE = "cvxr")
                                        }
                                       ,
                                        getI = function() {
                                                .Call('_cvxr_ProblemData__get_I', private$ptr, PACKAGE = "cvxr")
                                        }
                                       ,
                                        getJ = function() {
                                            .Call('_cvxr_ProblemData__get_J', private$ptr, PACKAGE = "cvxr")
                                        }
                                       ,
                                        getConstVec = function() {
                                            .Call('_cvxr_ProblemData__get_const_vec', private$ptr, PACKAGE = "cvxr")
                                        }
                                    ))
