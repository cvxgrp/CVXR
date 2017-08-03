## R6Solver class like solver.py in cvxpy
R6Solver <- R6::R6Class("R6Solver",
                        private = list(
                        ),
                        active = list(
                        ),
                        public = list(
                            name = function() {
                                ## The name of the solver
                            }
                           ,
                            import_solver = function() {
                                ## Imports the solver
                            }
                           ,
                            matrix_intf = function() {
                                ##The interface for matrices passed to the solver.
                            }
                           ,
                            vec_intf = function() {
                                ##The interface for vectors passed to the solver.
                            }
                           ,
                            split_constr = function(constr_map) {
                                ##Extracts the equality, inequality, and nonlinear constraints.
                                ##Parameters
                                ##----------
                                ##constr_map : dict
                                ##    A dict of the canonicalized constraints.
                            }
                           ,
                            choose_solver = function(constraints) {
                                ##Determines the appropriate solver.

                                ##Parameters
                                ##----------
                                ##constraints: list
                                ##    The list of canonicalized constraints.

                                ##Returns
                                ##-------
                                ##str
                                ##    The solver that will be used.
                                constr_map <- SymData.filter_constraints(constraints)
                                        # If no constraints, use ECOS.
                                if(length(constraints) == 0)
                                    return(ECOS())
                                        # If mixed integer constraints, use ECOS_BB.
                                else if(length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0)
                                    return(ECOS_BB())
                                        # If SDP, defaults to CVXOPT.
                                else if(length(constr_map[[SDP_MAP]]) > 0)
                                    return(CVXOPT())
                                        # Otherwise use ECOS
                                else
                                    return(ECOS())
                            }
                           ,

                            is_installed = function() {
                                ##Is the solver installed?
                                ##
                                tryCatch(requireNameSpace("xx", quietly = TRUE))
                                    self.import_solver()
                                return True
                                except ImportError:
                                           return False
                            }

                            toString = function() {
                                sparse <- self$sparse
                                if (sparse) {
                                    data <- paste(self$sparse_data, collapse=", ")
                                } else {
                                    data <- paste(self$dense_data, collapse=", ")
                                }
                                sprintf("LinOp(id=%s, type=%s, size=[%s], args=%s, sparse=%s, data=[%s])",
                                        self$get_id(),
                                        self$type,
                                        paste(self$size, collapse=", "),
                                        private$args$toString(),
                                        sparse,
                                        data)
                            }
                           ,
                            print = function() {
                                print(self$toString())
                            }

                        ))
