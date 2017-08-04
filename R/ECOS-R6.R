## ECOS class like ecos_int.py in cvxpy
ECOS_PACKAGE <- "ECOSolveR"

ECOS <- R6::R6Class("ECOS",
                    inherits = "R6Solver",
                    private = list(
                    ),
                    public = list(
                        initialize = function() {
                            private$name = "ECOS"
                            private$LP_CAPABLE <- TRUE
                            private$SOCP_CAPABLE <- TRUE
                            private$SDP_CAPABLE <- FALSE
                            private$EXPP_CAPABLE <- TRUE
                            private$MIP_CAPABLE <- FALSE
                        }
                        ,
                        import_solver = function() {
                            ## Imports the solver
                            requireNamespace(ECOS_PACKAGE, quietly = TRUE)
                        }
                       ,
                        matrix_intf = function() {
                            ##The interface for matrices passed to the solver.
                            DEFAULT_SPARSE_INTF
                        }
                       ,
                        vec_intf = function() {
                            ##The interface for vectors passed to the solver.
                            DEFAULT_INTF
                        }
                       ,
                        split_constr = function(constr_map) {
                            ##Extracts the equality, inequality, and nonlinear constraints.
                            ##Parameters
                            ##----------
                            ##constr_map : dict
                            ##    A dict of the canonicalized constraints.
                            list(eq_constr = constr_map[[EQ_MAP]],
                                 ineq_constr = constr_map[[LEQ_MAP]],
                                 nonlin_constr = list())
                        }
                       ,
                        solve = function(objective, constraints, cached_data, warm_start,
                                         verbose, solver_opts) {
                            self$import_solver()
                            data <- self$get_problem_data(objective, constraints, cached_data)
                            data[[DIMS]]['e'] <- data[[DIMS]][[EXP_DIM]]
                            results_dict <- ECOSolveR::ECOS_csolve(c = data[[C]],
                                                                   G = data[[G]],
                                                                   h = data[[H]],
                                                                   dims = data[[DIMS]],
                                                                   A = data[[A]],
                                                                   b = data[[B]],
                                                                   control = c(list(VERBOSE = verbose),
                                                                               solver_opts))
                            self$format_results(result_dict, data, cached_data)
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

                            new_results <- list()
                            status <- self$STATUS_MAP(results_dict["info"]["exitFlag"])
                            new_results[STATUS] <- status

                                        # Timing data
                            new_results[SOLVE_TIME] <- results_dict["info"]["timing"]["tsolve"]
                            new_results[SETUP_TIME] <- results_dict["info"]["timing"]["tsetup"]
                            new_results[NUM_ITERS] <- results_dict["info"]["iter"]

                            if(new_results[STATUS] %in% SOLUTION_PRESENT) {
                                primal_val <- results_dict['info']['pcost']
                                new_results[VALUE] <- primal_val + data[OFFSET]
                                new_results[PRIMAL] <- results_dict['x']
                                new_results[EQ_DUAL] <- results_dict['y']
                                new_results[INEQ_DUAL] <- results_dict['z']
                            }
                            new_results
                        }
                       ,
                        status_map = function(status) {
                            ## EXITCODES from ECOS
                            ## ECOS_OPTIMAL  (0)   Problem solved to optimality
                            ## ECOS_PINF     (1)   Found certificate of primal infeasibility
                            ## ECOS_DINF     (2)   Found certificate of dual infeasibility
                            ## ECOS_INACC_OFFSET (10)  Offset exitflag at inaccurate results
                            ## ECOS_MAXIT    (-1)  Maximum number of iterations reached
                            ## ECOS_NUMERICS (-2)  Search direction unreliable
                            ## ECOS_OUTCONE  (-3)  s or z got outside the cone, numerics?
                            ## ECOS_SIGINT   (-4)  solver interrupted by a signal/ctrl-c
                            ## ECOS_FATAL    (-7)  Unknown problem in solver

                            ## Map of ECOS status to CVXPY status.
                            if (status == 0) OPTIMAL
                            else if (status == 1) INFEASIBLE
                            else if (status == 2) UNBOUNDED
                            else if (status == 10) OPTIMAL_INACCURATE
                            else if (status == 11) INFEASIBLE_INACCURATE
                            else if (status == 12) UNBOUNDED_INACCURATE
                            else if (status %in% c(-1, -2, -3, -4, -7)) SOLVER_ERROR
                            else stop("ECOS status unrecognized: ", status)
                        })
                    ) ## End of R6Class

