## SCSS class like scs_int.py in cvxpy
SCS_PACKAGE <- "scs"

tri_to_full <- function(lower_tri, n, off_diag_scale = 1.0 / sqrt(2.0)) {
    ## Expands floor(n*(n+1)/2) lower triangular to full matrix, with off-diagonal entries
    ## scaled by scale defaulting to 1/sqrt(2)
    full <- matrix(0, nrow = n, ncol = n)
    for(col in 1:n) {
        for(row in col:n) {
            idx <- row - col + floor(n*(n+1)/2) - floor((n-col)*(n-col+1)/2)
            if(row != col) {
                full[row, col] <- lower_tri[idx] * off_diag_scale
                full[col, row] <- lower_tri[idx] * off_diag_scale
            } else
                full[row, col] <- lower_tri[idx]
        }
    }
    ##return(matrix(full, nrow = n^2))
    full
}


SCS <- R6::R6Class("SCS",
                   inherits = "R6Solver",
                   private = list(
                   ),
                   public = list(
                       initialize = function() {
                           private$name <- "SCS"
                           private$LP_CAPABLE <- TRUE
                           private$SOCP_CAPABLE <- TRUE
                           private$SDP_CAPABLE <- TRUE
                           private$EXPP_CAPABLE <- TRUE
                           private$MIP_CAPABLE <- FALSE
                       }
                      ,
                       import_solver = function() {
                           ## Imports the solver
                           requireNamespace(SCS_PACKAGE, quietly = TRUE)
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
                           list(eq_constr = c(constr_map[[EQ_MAP]],
                                              constr_map[[LEQ_MAP]]),
                                ineq_constr = list(),
                                nonlin_constr = list())
                       }
                      ,
                       solve = function(objective, constraints, cached_data, warm_start,
                                        verbose, solver_opts) {
                           self$import_solver()
                           data <- self$get_problem_data(objective, constraints, cached_data)
                           ## Set the options to be VERBOSE plus any user-specific options
                           solver_opts["verbose"] <- verbose
                           scs_args <- list(c = data[[C]], A = data[[A]], b = data[[B]])
                           ## If warm starting, add old primal and dual variables
                           solver_cache <- cached_data[self$name]
                           if(warm_start && !is.na(solver_cache@prev_result)) {
                               stop("Warm start currently unimplemented")
                               scs_args["x"] <- solver_cache@prev_result["x"]
                               scs_args["y"] <- solver_cache@prev_result["y"]
                               scs_args["s"] <- solver_cache@prev_result["s"]
                           }
                           ## results_dict <- do.call(scs::scs, c(scs_args, list(cone = data[[DIMS]]),
                           ##                         solver_opts))
                           results_dict <- scs::scs(A = data[[A]],
                                                    b = data[[B]],
                                                    obj = data[[C]],
                                                    cone = data[[DIMS]],
                                                    solver_opts)
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

                            solver_cache <- cached_data[self$name]
                            dims <- data[[DIMS]]
                            new_results <- list()
                            status <- self$status_map(results_dict["info"]["status"])
                            new_results[STATUS] <- status

                            ## Timing and iteration data
                            new_results[SOLVER_TIME] <- results_dict["info"]["solveTime"]/1000
                            new_results[SETUP_TIME] <- results_dict["info"]["setupTime"]/1000
                            new_results[NUM_ITERS] <- results_dict["info"]["iter"]

                            if (new_results[STATUS] %in% SOLUTION_PRESENT) {
                                ## Save previous result for possible future warm start
                                solver_cache@prev_result <- list(x = results_dict["x"],
                                                                 y = results_dict["y"],
                                                                 s = results_dict["s"])
                                primal_val <- results_dict["info"]["pobj"]
                                new_results[VALUE] <- primal_val + data[OFFSET]
                                new_results[PRIMAL] <- results_dict["x"]
                                new_results[EQ_DUAL] <- results_dict["y"][1:dims[EQ_DIM]]

                                y <- results_dict["y"][(dims[EQ_DIM]+1):length(results_dict["y"])]
                                old_sdp_sizes <- sum(sapply(dims[SDP_DIM], function(n) { floor(n*(n+1)/2) }))
                                new_sdp_sizes <- sum(dims[SDP_DIM]^2)
                                y_true <- rep(0, y@shape[1] + (new_sdp_sizes - old_sdp_sizes))
                                y_offset <- dims[LEQ_DIM] + sum(dims[SOC_DIM])
                                y_true_offset <- y_offset
                                y_true[1:y_true_offset] <- y[1:y_offset]

                                ## Expand SDP duals from lower triangular to full matrix, scaling off diagonal entries by 1/sqrt(2)
                                for (n in dims[SDP_DIM]) {
                                    tri <- y[y_offset:(y_offset + floor(n*(n+1)/2))]
                                    y_true[y_true_offset:(y_true_offset + n^2)] <- tri_to_full(tri, n)
                                    y_true_offset <- y_true_offset + n^2
                                    y_offset <- y_offset + floor(n*(n+1)/2)
                                }

                                y_true[(y_true_offset+1):length(y_true)] <-y[(y_offset+1):length(y)]
                                new_results[INEQ_DUAL] <- y_true
                            } else {
                                ## No result to save
                                solver_cache@prev_result <- NA
                            }
                            new_results
                        }
                       ,
                        status_map = function(status) {
                            ## EXITCODES from SCS
                            if (status == "Solved") OPTIMAL
                            else if (status == "Solved/Inaccurate") OPTIMAL_INACCURATE
                            else if (status == "Unbounded") UNBOUNDED
                            else if (status == "Unbounded/Inaccurate") UNBOUNDED_INACCURATE
                            else if (status == "Infeasible") INFEASIBLE
                            else if (status == "Infeasible/Inaccurate") INFEASIBLE_INACCURATE
                            else if (status %in% c("Failure", "Indeterminate")) SOLVER_ERROR
                            else stop("SCS status unrecognized: ", status)
                        })
                    ) ## End of R6Class

