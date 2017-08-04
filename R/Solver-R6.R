## R6Solver class like solver.py in cvxpy
R6Solver <- R6::R6Class("R6Solver",
                        private = list(
                            name = NA,
                            LP_CAPABLE = NA,
                            SOCP_CAPABLE = NA,
                            SDP_CAPABLE = NA,
                            EXP_CAPABLE = NA,
                            MIP_CAPABLE = NA,
                            reject_problem = function(reason) {
                                ##Raise an error indicating that the solver cannot solve a problem.
                                ##Parameters
                                ##----------
                                ## reason : str
                                ##A short description of the reason the problem cannot be solved by
                                ##this solver.
                                ##Raises
                                ##------
                                ##cvxpy.SolverError
                                ##An error explaining why the problem could not be solved.
                                message = sprintf("The solver %s cannot solve the problem because %s.",
                                                  self$name, reason)
                                stop(message)
                            }
                           ,
                            noncvx_id_to_idx = function(dims, var_offsets, var_sizes) {
                                bool_idx <- lapply(dims[BOOL_IDS], function(var_id) {
                                    offset <- var_offsets[var_id]
                                    size <- var_sizes[var_id]
                                    offset + seq(1, size[1]*size[2], by = 1)
                                })

                                int_idx <- lapply(dims[INT_IDS], function(var_id) {
                                    offset <- var_offsets[var_id]
                                    size <- var_sizes[var_id]
                                    offset + seq(1, size[1]*size[2], by = 1)
                                })

                                list(bool_idx = bool_idx, int_idx = int_idx)
                            }
                        ),
                        active = list(
                            name = function(value) {
                                if (missing(value)) {
                                    private$name <- value
                                } else {
                                    private$name
                                }
                            }
                           ,
                            LP_CAPABLE = function(value) {
                                if (missing(value)) {
                                    private$LP_CAPABLE <- value
                                } else {
                                    private$LP_CAPABLE
                                }
                            }
                           ,
                            SOCP_CAPABLE = function(value) {
                                if (missing(value)) {
                                    private$SOCP_CAPABLE <- value
                                } else {
                                    private$SOCP_CAPABLE
                                }
                            }
                           ,
                            SDP_CAPABLE = function(value) {
                                if (missing(value)) {
                                    private$SDP_CAPABLE <- value
                                } else {
                                    private$SDP_CAPABLE
                                }
                            }
                           ,
                            EXP_CAPABLE = function(value) {
                                if (missing(value)) {
                                    private$EXP_CAPABLE <- value
                                } else {
                                    private$EXP_CAPABLE
                                }
                            }
                           ,
                            MIP_CAPABLE = function(value) {
                                if (missing(value)) {
                                    private$MIP_CAPABLE <- value
                                } else {
                                    private$MIP_CAPABLE
                                }
                            }
                        ),
                        public = list(
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
                                ## If no constraints, use ECOS.
                                if (length(constraints) == 0) {
                                    ECOS$new()
                                    ## If mixed integer constraints, use ECOS_BB.
                                } else if (length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0) {
                                    ECOS_BB$new()
                                    ## If SDP, defaults to CVXOPT.
                                } else if (length(constr_map[[SDP_MAP]]) > 0) {
                                    SCS$new()
                                    ## Otherwise use ECOS
                                } else {
                                    ECOS$new()
                                }
                            }
                           ,
                            is_installed = function() {
                                ##Is the solver installed?
                                ##
                                tryCatch(self$import_solver(),
                                         error = function(e) FALSE)
                            }
                           ,
                            validate_solver = function(constraints) {
                                ##Raises an exception if the solver cannot solve the problem.
                                ##Parameters
                                ##----------
                                ##constraints: list
                                ##    The list of canonicalized constraints.
                                ## Check the solver is installed.
                                if (!self$is_installed())
                                    stop(sprintf("The solver %s is not installed.", self$name()))
                                ## Check the solver can solve the problem.
                                constr_map <- SymData.filter_constraints(constraints)

                                if ( (length(constr_map[[BOOL_MAP]]) > 0 ||
                                      length(constr_map[[INT_MAP]]) > 0) && !self$MIP_CAPABLE) {
                                    private$reject_problem("it cannot solve mixed-integer problems")
                                } else if (length(constr_map[[SDP_MAP]]) > 0 && !self$SDP_CAPABLE) {
                                    private$reject_problem("it cannot solve semidefinite problems")
                                } else if (length(constr_map[[EXP_MAP]]) > 0 && !self$EXP_CAPABLE) {
                                    private$reject_problem("it cannot solve exponential cone problems")
                                } else if (length(constr_map[[SOC_MAP]]) > 0 && !self$SOCP_CAPABLE) {
                                    private$reject_problem("it cannot solve second-order cone problems")
                                } else if (length(constraints) == 0 && self$name() %in% c("SCS", "GLPK")) {
                                    private$reject_problem("it cannot solve unconstrained problems")
                                } else {
                                    invisible(TRUE)
                                }
                            }
                           ,
                            validate_cache = function(objective, constraints, cached_data) {
                                ##Clears the cache if the objective or constraints changed.
                                ##Parameters
                                ##----------
                                ##objective : LinOp
                                ##The canonicalized objective.
                                ##constraints : list
                                ##The list of canonicalized cosntraints.
                                ##cached_data : dict
                                ##A map of solver name to cached problem data.
                                prob_data <- cached_data[[self$name]]
                                if (!is.null(prob_data@sym_data) &&
                                    (!isTRUE(all.equal(objective, prob_data@sym_data@objective)) ||
                                     !isTRUE(all.equal(constraints, prob_data@sym_data@constraints)))) {
                                    prob_data@sym_data <- NULL
                                    prob_data@matrix_data <- NULL
                                }
                                cached_data[[self$name]] <- prob_data
                                cached_data
                            }
                           ,
                            get_sym_data = function(objective, constraints, cached_data) {
                                ##Returns the symbolic data for the problem.
                                ##Parameters
                                ##----------
                                ##objective : LinOp
                                ##    The canonicalized objective.
                                ##constraints : list
                                ##    The list of canonicalized constraints.
                                ##cached_data : dict
                                ##    A map of solver name to cached problem data.
                                ##Returns
                                ##-------
                                ##SymData
                                ##    The symbolic data for the problem.
                                cached_data <- self$validate_cache(objective, constraints, cached_data)
                                prob_data <- cached_data[[self$name]]
                                if (is.null(prob_data@sym_data))
                                    prob_data@sym_data <- SymData(objective, constraints, solver)
                                cached_data[[self$name]] <- prob_data
                                cached_data
                            }
                           ,
                            get_problem_data = function(objective, constraints, cached_data) {
                                ##Returns the numeric data for the problem.
                                ##Parameters
                                ##----------
                                ##objective : LinOp
                                ##    The canonicalized objective.
                                ##constraints : list
                                ##    The list of canonicalized cosntraints.
                                ##cached_data : dict
                                ##    A map of solver name to cached problem data.
                                ##Returns
                                ##-------
                                ##SymData
                                ##    The symbolic data for the problem.
                                cached_data <- self$get_sym_data(objective, constraints, cached_data)
                                sym_data <- cached_data[[self$name]]@sym_data
                                cached_data <- self$get_matrix_data(objective, constraints, cached_data)
                                matrix_data <- cached_data[[self$name]]@matrix_data

                                data <- list()
                                obj <- get_objective(matrix_data)
                                eq <- get_eq_constr(matrix_data)
                                ineq <- get_ineq_constr(matrix_data)

                                data[[C]] <- obj[[1]]
                                data[[OFFSET]] <- obj[[2]]
                                data[[A]] <- eq[[1]]
                                data[[B]] <- eq[[2]]
                                data[[G]] <- ineq[[1]]
                                data[[H]] <- ineq[[2]]
                                data[[DIMS]] <- sym_data@dims

                                conv_idx <- private$noncvx_id_to_idx(data[[DIMS]],
                                                                     sym_data@var_offsets,
                                                                     sym_data@var_sizes)
                                data[[BOOL_IDX]] <- conv_idx$bool_idx
                                data[[INT_IDX]] <- conv_idx$int_idx
                                data
                            }
                           ,
                            nonlin_constr = function() {
                                ##Returns whether nonlinear constraints are needed.
                                False
                            }
                           ,
                            is_mip = function(data) {
                                ##Is the problem a mixed integer program?
                                length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
                            }
                           ,
                            format_results = function(results_dict, data, cached_data) {
                                ##Converts the solver output into standard form.
                                ##Parameters
                                ##----------
                                ##results_dict : dict
                                ##    The solver output.
                                ##data : dict
                                ##    Information about the problem.
                                ##cached_data : dict
                                ##    A map of solver name to cached problem data.
                                ##Returns
                                ##-------
                                ##dict
                                ##    The solver output in standard form.
                            }
                        ))

