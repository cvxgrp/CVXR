## CVXcanon class shadowing CVXcannon.cpp code, exposing merely the two
## build_matrix methods
CVXcanon <- R6::R6Class("CVXcanon",
                        private = list(
                            ptr = NA
                        )
                       ,
                        active = list(
                        )
                       ,
                        public = list(
                            initialize = function() {
                            }
                           ,
                            getXPtr = function() {
                                private$ptr
                            }
                           ,
                            build_matrix = function(constraints, id_to_col, constr_offsets) {
                                ## constraints is a vector of Linops (LinOpVector-R6 in R)
                                ## id_to_col is an integer vector with names that are
                                ## integers converted to chacracters
                                ## constr_offsets is a standard integer vector in R
                                ## cat("Linvec\n")
                                ## constraints$print()
                                ## cat("id_to_col\n")
                                ## print(id_to_col)

                                if (missing(constr_offsets)) {
                                    objPtr <- .Call('_CVXR_build_matrix_0',
                                                    constraints$getXPtr(),
                                                    id_to_col,
                                                    PACKAGE = 'CVXR')
                                } else {
                                    objPtr <- .Call('_CVXR_build_matrix_1',
                                                    constraints$getXPtr(),
                                                    id_to_col,
                                                    constr_offsets,
                                                    PACKAGE = 'CVXR')
                                }
                                ##cat("Instantiating ProblemData-R6", "\n")
                                ##browser()
                                CVXcanon.ProblemData$new(objPtr)
                            }
                        ))
