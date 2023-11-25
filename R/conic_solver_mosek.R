#'
#' Given a problem returns a PSD constrain
#'
#' @param problem A \linkS4class{Problem} object.
#' @param c A vector of coefficients.
#' @return Returns an array G and vector h such that the given constraint is
#' equivalent to G*z <=_{PSD} h.
psd_coeff_offset <- function(problem, c) {
  extractor <- CoeffExtractor(InverseData(problem))
  tmp <- affine(extractor, expr(c))
  A_vec <- tmp[[1]]
  b_vec <- tmp[[2]]
  G <- -A_vec
  h <- b_vec
  dim <- nrow(expr(c))
  return(list(G, h, dim))
}

#' An interface for the MOSEK solver.
#'
#' @name MOSEK-class
#' @aliases MOSEK
#' @rdname MOSEK-class
#' @export
MOSEK <- setClass("MOSEK",
                  slots = list(DUAL_EXP_CONE_ORDER = "integer",
                               MI_SUPPORTED_CONSTRAINTS = "character"),
                  prototype = list(EXP_CONE_ORDER = c(2L, 1L, 0L),
                                   DUAL_EXP_CONE_ORDER = c(0L, 1L, 2L),
                                   MIP_CAPABLE = TRUE,
                                   SUPPORTED_CONSTRAINTS = c(ConicSolver()@SUPPORTED_CONSTRAINTS, "SOC", "PSDConstraint", "ExpCone", "PowCone3D"),
                                   MI_SUPPORTED_CONSTRAINTS = c(ConicSolver()@SUPPORTED_CONSTRAINTS, "SOC", "ExpCone", "PowCone3D")),
                  contains = "ConicSolver")

#' @param object,x A \linkS4class{MOSEK} object.
#' @describeIn MOSEK Imports the solver.
#' @importFrom utils packageDescription
setMethod("import_solver", "MOSEK", function(solver) {
  if (requireNamespace("Rmosek", quietly = TRUE)) {
    if (is.null(utils::packageDescription("Rmosek")$Configured.MSK_VERSION)) {
      message("Rmosek package installed, but not configured per instructions therein!")
      FALSE
    } else {
      TRUE
    }
    ## TODO: Add exponential cone support.
  } else {
    message("Rmosek not installed.")
    FALSE
  }
})

#' @describeIn MOSEK Returns the name of the solver.
setMethod("name", "MOSEK", function(x) { MOSEK_NAME })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn MOSEK Can MOSEK solve the problem?
setMethod("accepts", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  import_solver(object)
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  sp <- object@SUPPORTED_CONSTRAINTS
  for(constr in problem@constraints) {
    if(!inherits(constr, sp))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

#' Construct a linear operator to multiply by PSD constraint coefficients.
#' @details
#' Special cases PSD constraints, as MOSEK expects constraints to be
#' imposed on solely the lower triangular part of the variable matrix.
#' This function differs from the SCS version only in the fact that it does not apply the `\sqrt(2)` scaling on off-diagonal entries. This difference is necessary based on how we implement the `MOSEK$bar_data`.
#' @param constr a list of constraints
#' @return a linear operator to multiply by PSD constraint coefficients.
#' @importFrom Matrix sparseMatrix
setMethod("psd_format_mat", "MOSEK",  function(solver, constr) {
  rows <- cols <- nrow(constr)
  ind_seq <- seq_len(rows)
  entries <- rows * (cols + 1) %/% 2
  val_array <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                    dims = c(rows, cols))
  val_tri <- lower.tri(val_array, diag = TRUE)
  col_arr <- which(val_tri)
  val_array[val_tri] <- 1.0
  diag(val_array) <- 1.0
  n <- rows * cols
  scaled_lower_tri <- Matrix::sparseMatrix(i = seq_len(entries), j = col_arr,
                                           x = val_array@x, dims = c(entries, n))
  idx <- seq_len(n)
  val_symm <- rep(0.5 , 2 * n)
  tmp <- ind_seq - 1L
  res <- rows * tmp
  row_symm <- c(idx, idx)
  col_symm <- c(idx, unlist(lapply(ind_seq, `+`, res)))
  symm_matrix <- Matrix::sparseMatrix(i = row_symm, j = col_symm, x = val_symm,
                                      dims = c(n, n))
  scaled_lower_tri %*% symm_matrix
})


MOSEK.bar_data <- function(A_psd, c_psd, K) {
  # TODO: investigate how to transform or represent "A_psd" so that the following
  #  indexing and slicing operations are computationally cheap. Or just rewrite
  #  the function so that explicit slicing is not needed on a SciPy sparse matrix.
  n <- nrow(A_psd)
  iseq <- seq_len(n)
  c_bar_data <- A_bar_data <- list()
  idx <- 0
  psd_variables <- K[[Cone2Cone.PSD]]
  for (j in seq_along(psd_variables)) {
    dim <- psd_variables[j]  # psd variable index j.
    vec_len <- dim * (dim + 1) / 2
    A_block <- A_psd[ , idx:(idx + vec_len - 1)]  ## is the -1 correct?
    # ^ each row specifies a linear operator on PSD variable.
    for (i in iseq) {
      A_row <- A_block[i, ]
      if (nnzero(A) >= 0) {
        A_row <- as.numeric(A_row)
        # A_row defines a symmetric matrix by where the first "order" entries
        #   gives the matrix's first column, the second "order-1" entries gives
        #   the matrix's second column (diagonal and blow), and so on.
        vlt <- vectorized_lower_tri_to_triples(A_row, dim)
        # TODO: replace the above function with something that only reads the nonzero
        #  entries. I.e. return *actual* sparse matrix data, rather than data for a dense
        #  matrix stated in a sparse format.
        A_bar_data <- append(A_bar_data, list(i = i, j = j, rows = vlt$rows, cols = vlt$cols, vals = vlt$vals))
      }
    }
    c_block <- c_psd[idx:idx + vec_len]
    vlt <- vectorized_lower_tri_to_triples(c_block, dim)
    c_bar_data <- append(c_bar_data, list(j = j, rows = vlt$rows, cols = vlt$cols, vals = vlt$vals))
    idx <- idx + vec_len
  }
  list(A = A_bar_data, c = c_bar_data)
}

#' @describeIn MOSEK Returns a new problem and data for inverting the new solution.
setMethod("perform", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
    if (! problem@formatted) {
        problem <- format_constr(object, problem, object@EXP_CONE_ORDER)
    }
    if (length(problem@x@boolean_idx) > 0 || length(problem@x@integer_idx) > 0) {
        res <- Slacks.perform(problem, list(Cone2Cone.NONNEG))
        data <- res$data
        inv_data <- res$inv_data
    } else {
        res <- Dualize.perform(problem)
        data <- res$data
        inv_data <- res$inv_data
        ## need to do more to handle SDP.
        A <- data[[A_KEY]]; c <- data[[C_KEY]]; K <- data$K_dir;
        num_psd <- length(K[[Cone2Cone.PSD]])
        if (num_psd > 0) {
            idx <- K[[Cone2Cone.FREE]] + K[[Cone2Cone.NONNEG]] + sum(K[[Cone2Cone.SOC]])
            total_psd <- sum(unlist(lapply(K[[Cone2Cone.PSD]], function(d) d * (d + 1) / 2)))
            ind_seq <- seq.int(from = idx, length.out = total_psd)
            A_psd <- A[ , ind_seq ]
            c_psd <- c[ind_seq]
            iseq <- seq_len(idx)
            if (K[[Cone2Cone.DUAL_EXP]] == 0) {
                data[[A_KEY]] = A[ , iseq]
                data[[C_key]] = c[iseq]
            } else {
                iseq_tail <- seq.int(from = idx + total_psd, to = ncol(A))
                data[[A_KEY]] <- cbind2(A[, iseq ], A[, iseq_tail])
                data[[C_KEY]] <- c[c(iseq, iseq_tail)]
            }
            bar_data <- MOSEK.bar_data(A_psd, c_psd, K)
            data[['A_bar_data']] = bar_data$A
            data[['c_bar_data']] = bar_data$c
        } else {
            data[['A_bar_data']] = list()
            data[['c_bar_data']] = list()
        }
    }
    data[[PARAM_PROB]] <- problem
    list(data = data, inv_data = inv_data)
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn MOSEK Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "MOSEK", function(object, data, warm_start,
                                              verbose, feastol, reltol, abstol,
                                              num_iter,
                                              solver_opts, solver_cache) {

    if (missing(solver_cache)) solver_cache <- new.env(parent = emptyenv())
  ## Check if the CVXR standard form has zero variables. If so,
  ## return a trivial solution. This is necessary because MOSEK
  ## will crash if handed a problem with zero variables.

  c <- data[[C_KEY]]
  if (data$dualized) {
    if (length(data[[C_KEY]]) == 0 && length(data[['c_bar_data']]) == 0) {
      # primal problem was unconstrained minimization of a linear function.
      if (norm(data[[B_KEY]], "2") > 0) {
        return( Solution(INFEASIBLE, -Inf, NA, NA, list()) )
      } else {
        return( Solution(OPTIMAL, 0.0, list(), list(EQ_DUAL = data[[B_KEY]]), list()) )
      }
    } else {
      task <- MOSEK.build_dualized_task(data)
    }
  } else {
    if (length(data[[C_KEY]]) == 0) {
      return( Solution(OPTIMAL, 0.0, list(), list(), list()) )
    } else {
      task <- MOSEK.build_slack_task(data)
    }
  }

  # Set parameters, optimize the Mosek Task, and return the result.
  ## save_file = MOSEK.handle_options(env, task, verbose, solver_opts)
  ## if save_file:
  ##      task.writedata(save_file)

  ## TODO:  Modify what's below to do warm solve and also take in options specified! See old code.
  result <- mosek(task);

  ## r <- Rmosek::mosek(prob, list(verbose = verbose, usesol = solver_opts$usesol,
  ##                               useparam = solver_opts$useparam, soldetail = solver_opts$soldetail,
  ##                               getinfo = solver_opts$getinfo, writebefore = solver_opts$writebefore,
  ##                               writeafter = solver_opts$writeafter))


  ## if (verbose) {
  ##   task.solutionsummary(mosek.streamtype.msg)
  ## }
  return(result)

})

MOSEK.build_dualized_task <- function(data)  {
  ## This function assumes "data" is formatted according to MOSEK.apply when the problem
  ## features no integer constraints. This dictionary should contain keys s.C, s.A, s.B,
  ## 'K_dir', 'c_bar_data' and 'A_bar_data'.

  ## If the problem has no PSD constraints, then we construct a Task representing

  ##    max{ c.T @ x : A @ x == b, x in K_dir }

  ## If the problem has PSD constraints, then the Task looks like

  ##    max{ c.T @ x + c_bar(X_bars) : A @ x + A_bar(X_bars) == b, x in K_dir, X_bars PSD }

  ## In the above formulation, c_bar is effectively specified by a list of appropriately
  ## formatted symmetric matrices (one symmetric matrix for each PSD variable). A_bar
  ## is specified a collection of symmetric matrix data indexed by (i, j) where the j-th
  ## PSD variable contributes a certain scalar to the i-th linear equation in the system
  ## "A @ x + A_bar(X_bars) == b".

  # problem data
  c <- data[[C_KEY]]; A <- data[[A_KEY]]; b <- data[[B_KEY]]; K <- data[['K_dir']];
  n <- nrow(A); m <- ncol(A);

  ## No task for Rmosek, but declares the number of variables in the model.
  task <- list()

  ## Need to expand prob$c as well to match this dimension
  ## task.appendvars(m)
  ## o <- numeric(m)
  ## task.putvarboundlist(np.arange(m, dtype=int), [mosek.boundkey.fr] * m, o, o)
  ## task.appendcons(n)
  blx <- bux <- rep(0, m)
  ## Remember to create, later below,
  ## task$bx <- rbind(blx = blx, bux = bux)

  ## objective
  ## task.putclist(np.arange(c.size, dtype=int), c)
  ## task.putobjsense(mosek.objsense.maximize)
  task$c <- c
  task$sense <- "max"

  ## equality constraints
  ## rows, cols, vals = sp.sparse.find(A)
  ## task.putaijlist(rows.tolist(), cols.tolist(), vals.tolist())
  task$A <- A

  ## task.putconboundlist(np.arange(n, dtype=int), [mosek.boundkey.fx] * n, b, b)
  task$bc <- rbind(blx = b, bux = b)

  ## conic constraints
  cones <- list()
  idx <- K[[Cone2Cone.FREE]]
  num_pos <- K[[Cone2Cone.NONNEG]]
  if (num_pos > 0) {
    ## o = np.zeros(num_pos)
    ## task.putvarboundlist(np.arange(idx, idx + num_pos, dtype=int),
    ##                      [mosek.boundkey.lo] * num_pos, o, o)
    iseq <- seq.int(from = idx + 1, length.out = num_pos)
    blx[iseq] <- bux[iseq] <- 0
    idx <- idx + num_pos
  }

  ## We can now finalize variable bounds
  task$bx <- rbind(blx = blx, bux = bux)


  cone_sizes <- K[[Cone2Cone.SOC]]
  num_soc <- length(cone_sizes)
  if (num_soc > 0) {
    ## Accumulate cones
    cones <- c(cones,
               lapply(one_sizes, function(dim) list("QUAD", dim, NULL)))
    idx <- idx + sum(as.integer(cone_sizes))
    ## cones = [mosek.conetype.quad] * num_soc
    ## task.appendconesseq(cones, [0] * num_soc, K[Cone2Cone.SOC], idx)
    ## idx += sum(K[Cone2Cone.SOC])
  }

  num_dexp <- K[[Cone2Cone.DUAL_EXP]]
  if (num_dexp > 0) {
    cones <- c(cones,
               lapply(cone_sizes, function(dummy) list("DEXP", 3, NULL)))
    idx <- idx + 3 * num_dexp

    ## cones = [mosek.conetype.dexp] * num_dexp
    ## task.appendconesseq(cones, [0] * num_dexp, [3] * num_dexp, idx)
    ## idx += 3 * num_dexp
  }

  cone_sizes <- K[[Cone2Cone.DUAL_POW3D]]
  num_dpow <- length(cone_sizes)
  if (num_dpow > 0) {
    cones <- c(cones,
               lapply(cone_sizes, function(dummy) list("DPOW", 3, NULL)))
    idx <- idx + 3 * num_dpow

    ## cones = [mosek.conetype.dpow] * num_dpow
    ## task.appendconesseq(cones, K[Cone2Cone.DUAL_POW3D], [3] * num_dpow, idx)
    ## idx += 3 * num_dpow
  }

  ## We can now finalize the cone specification matrix
  task$cones <- do.call(matrix, c(cones, list(nrow = 3)))

  psd_dims <- K[[Cone2Cone.PSD]]
  num_psd <- length(psd_dims)
  if (num_psd > 0) {
    ## task.appendbarvars(K[Cone2Cone.PSD])
    ## psd_dims <- K[[Cone2Cone.PSD]]
    i_vec <- j_vec <- k_vec <- v_vec <- c()
    for (item in data['A_bar_data']) {
      n  <- length(item$vals)
      i_vec  <- c(i_vec, rep(item$i, n))
      j_vec  <- c(j_vec, rep(item$j, n))
      k_vec  <- c(k_vec, item$rows)
      l_vec  <- c(l_vec, item$cols)
      v_vec  <- c(v_vec, item$vals)
    }
    task$barA$i  <- i_vec; task$barA$j  <- j_vec;
    task$barA$k  <- k_vec; task$barA$l  <- l_vec;
    task$barA$v  <- v_vec;

    j_vec <- k_vec <- v_vec <- c()
    for (item in data['c_bar_data']) {
      n  <- length(item$vals)
      j_vec  <- c(j_vec, rep(item$j, n))
      k_vec  <- c(k_vec, item$rows)
      l_vec  <- c(l_vec, item$cols)
      v_vec  <- c(v_vec, item$vals)
    }
    task$barc$j  <- j_vec; task$barc$k <- k_vec;
    task$barc$l  <- l_vec; task$barc$v  <- v_vec;
  }
  return(task)
}

MOSEK.build_slack_task <- function(data) {
  ## This function assumes "data" is formatted by MOSEK.apply, and is only intended when
  ## the problem has integer constraints. As of MOSEK version 9.2, MOSEK does not support
  ## mixed-integer SDP. This implementation relies on that fact.

  ## "data" is a dict, keyed by s.C, s.A, s.B, 'K_dir', 'K_aff', s.BOOL_IDX and s.INT_IDX.
  ## The data 'K_aff' corresponds to constraints which MOSEK accepts as "A @ x <=_{K_aff}"
  ## (so-called "affine"  conic constraints), in contrast with constraints which must be stated
  ## as "x in K_dir" ("direct" conic constraints). As of MOSEK 9.2, the only allowed K_aff is
  ## the zero cone and the nonnegative orthant. All other constraints must be specified in a
  ## "direct" sense.

  ## The returned Task represents

  ##     min{ c.T @ x : A @ x <=_{K_aff} b,  x in K_dir, x[bools] in {0,1}, x[ints] in Z }.

  K_aff <- data[['K_aff']]
  # K_aff keyed by a2d.ZERO, a2d.NONNEG
  c  <- data[[C_KEY]]; A <- data[[A_KEY]];  b <- data[[B_KEY]];
  # The rows of (A, b) go by a2d.ZERO and then a2d.NONNEG
  K_dir <- data['K_dir']
  # Components of the vector "x" are constrained in the order
  # a2d.FREE, then a2d.SOC, then a2d.EXP. PSD is not supported.
  d <- dim(A); m  <- d[1]; n  <- d[2];
  ## No task for Rmosek, but declares the number of variables in the model.
  task <- list()
  # task.appendvars(n)
  # o = np.zeros(n)
  # task.putvarboundlist(np.arange(n, dtype=int), [mosek.boundkey.fr] * n, o, o)
  #task.appendcons(m)
  blx  <- bux  <- rep(0, n)
  ## Remember to create, later below,
  ## task$bx <- rbind(blx = blx, bux = bux)

  ## objective
  ## task.putclist(np.arange(n, dtype=int), c)
  ## task.putobjsense(mosek.objsense.minimize)
  task$c  <- c
  task$sense  <- "min"
  ## elementwise constraints
  ## rows, cols, vals = sp.sparse.find(A)
  ## task.putaijlist(rows, cols, vals)
  task$A  <- A
  ## eq_keys = [mosek.boundkey.fx] * K_aff[a2d.ZERO]
  ## ineq_keys = [mosek.boundkey.up] * K_aff[a2d.NONNEG]
  ## task.putconboundlist(np.arange(m, dtype=int), eq_keys + ineq_keys, b, b)
  task$bc <- rbind(blx = b, bux = b)

  ## conic constraints
  cones <- list()
  idx <- K[[Cone2Cone.FREE]]
  cone_sizes  <- K_dir[[Cone2Cone.SOC]]
  num_soc <- length(cone_sizes)
  if (num_soc > 0) {
    cones <- c(cones,
               lapply(cone_sizes, function(dim) list("QUAD", dim, NULL)))
    idx <- idx + sum(as.integer(cone_sizes))
    ## conetypes = [mosek.conetype.quad] * num_soc
    ## task.appendconesseq(conetypes, [0] * num_soc, K_dir[a2d.SOC], idx)
    ## idx += sum(K_dir[a2d.SOC])
  }

  cone_sizes <- K_dir[[Cone2Cone.EXP]]
  num_exp  <- length(cone_sizes)
  if (num_exp > 0) {
    cones <- c(cones,
               lapply(cone_sizes, function(dummy) list("PEXP", 3, NULL)))
    idx <- idx + 3 * num_exp

    ## conetypes = [mosek.conetype.pexp] * num_exp
    ## task.appendconesseq(conetypes, [0] * num_exp, [3] * num_exp, idx)
    ## idx += 3*num_exp
  }

  cone_sizes <- K_dir[[Cone2Cone.POW3D]]
  num_pow  <- length(cone_sizes)
  if (num_pow > 0) {
    cones <- c(cones,
               lapply(cone_sizes, function(dummy) list("DPOW", 3, NULL)))
    idx <- idx + 3 * num_pow
    ## conetypes = [mosek.conetype.ppow] * num_pow
    ## task.appendconesseq(conetypes, K_dir[a2d.POW3D], [3] * num_pow, idx)
    ## idx += 3*num_pow
  }

  ## integrality constraints
  bools  <- data[[BOOL_IDX]]
  num_bool <- length(bools)
  ints <- data[[INT_IDX]]
  num_int  <- length(ints)
  task$intsub  <- c(bools, ints)

  ## We can now finalize variable bounds
  blx[bools]  <- 0; bux[bools] <- 1;
  task$bx <- rbind(blx = blx, bux = bux)
  ## vartypes = [mosek.variabletype.type_int] * (num_bool + num_int)
  ## task.putvartypelist(data[s.INT_IDX] + data[s.BOOL_IDX], vartypes)
  ## if num_bool > 0:
  ##      task.putvarboundlist(data[s.BOOL_IDX], [mosek.boundkey.ra] * num_bool,
  ##                           [0] * num_bool, [1] * num_bool)
  return(task)
}

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn MOSEK Returns the solution to the original problem given the inverse_data.
setMethod("invert", "MOSEK", function(object, solution, inverse_data) {
  ## REMOVE LATER
  ## results  <- solution
  ##    has_attr <- !is.null(mosek.solsta$near_optimal)
  ## We ignore MOSEK 8.1 and below.
  status_map <- function(status) {
      status  <- tolower(status)
      if(status %in% c("optimal", "integer_optimal"))
          return(OPTIMAL)
      ##        else if(status %in% c("prim_feas", "near_optimal", "near_integer_optimal"))
      ##            return(OPTIMAL_INACCURATE)
      else if(status == "prim_infeas_cer" || status == "primal_infeasible_cer") { #Documentation says it's this, but docs also say it spits out dual_infeas_cer, which is wrong
        #check later
          if(!is.null(attributes(status))) #check if status has any attributes, hasattr in python
              return(INFEASIBLE)
          else
              return(INFEASIBLE)
      } else if(status == "dual_infeasible_cer") {
          if(!is.null(attributes(status)))
              return(UNBOUNDED_INACCURATE)
          else
              return(UNBOUNDED)
      }  else
          return(SOLVER_ERROR)
  }

  ##env <- results$env
  ##task <- results$task
  ## Naras: FIX solver_opts
  solver_opts <- solution$solver_options

  if(inverse_data$integer_variables)
      sol <- solution$sol$int
  else if(!is.null(solver_opts$bfs) && solver_opts$bfs && inverse_data$is_LP)
      sol <- solution$sol$bas   # The basic feasible solution.
  else
      sol <- solution$sol$itr   # The solution found via interior point method.

  problem_status <- sol$prosta
  solution_status <- sol$solsta

  if(is.na(solution$response$code))
    status <- SOLVER_ERROR
  else
    status <- status_map(solution_status)

  ## For integer problems, problem status determines infeasibility (no solution).
  ##  if(sol == mosek.soltype.itg && problem_status == mosek.prosta.prim_infeas)
  ## Using reference https://docs.mosek.com/9.0/rmosek/accessing-solution.html
  if(inverse_data$integer_variables && (problem_status == "MSK_PRO_STA_PRIM_INFEAS" || problem_status == "PRIMAL_INFEASIBLE"))
      status <- INFEASIBLE

  if(status %in% SOLUTION_PRESENT) {
                                      # Get objective value.
      opt_val <- sol$pobjval + inverse_data[[OBJ_OFFSET]]
                                      # Recover the CVXR standard form primal variable.
      ## z <- rep(0, inverse_data$n0)
      ## task.getxxslice(sol, 0, length(z), z)
      primal_vars <- list()
      primal_vars[[as.character(inverse_data[[object@var_id]])]] <- sol$xx
      ## Recover the CVXR standard form dual variables.
      ## if(sol == mosek.soltype.itn)
      if (inverse_data$integer_variables) {
        dual_var_ids <- sapply(c(inverse_data$suc_slacks, inverse_data$y_slacks, inverse_data$snx_slacks, inverse_data$psd_dims), function(slack) { slack[[1L]] })
        dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
        names(dual_vars) <- dual_var_ids
      } else
        dual_vars <- MOSEK.recover_dual_variables(sol, inverse_data)

  } else {
      if(status == INFEASIBLE)
          opt_val <- Inf
      else if(status == UNBOUNDED)
          opt_val <- -Inf
      else
          opt_val <- NA_real_
      vid <-
      primal_vars <- list()
      primal_vars[[as.character(inverse_data[[object@var_id]])]] <- NA_real_
      dual_var_ids <- sapply(c(inverse_data$suc_slacks, inverse_data$y_slacks, inverse_data$snx_slacks, inverse_data$psd_dims), function(slack) { slack[[1L]] })
      dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
      names(dual_vars) <- dual_var_ids
  }

  ## Store computation time.
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$dinfo$OPTIMIZER_TIME

  ## Delete the MOSEK Task and Environment
  ##task.__exit__(NA, NA, NA)
  ##env.__exit__(NA, NA, NA)

  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

#'
#' Recovers MOSEK solutions dual variables
#'
#' @param sol List of the solutions returned by the MOSEK solver.
#' @param inverse_data A list of the data returned by the perform function.
#' @return A list containing the mapping of CVXR's \linkS4class{Constraint}
#' object's id to its corresponding dual variables in the current solution.
MOSEK.recover_dual_variables <- function(sol, inverse_data) {
  dual_vars <- list()

  ## Dual variables for the inequality constraints.
  suc_len <- ifelse(length(inverse_data$suc_slacks) == 0, 0, sum(sapply(inverse_data$suc_slacks, function(val) { val[[2]] })))
  if(suc_len > 0) {
      ## suc <- rep(0, suc_len)
      ## task.getsucslice(sol, 0, suc_len, suc)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$suc[seq_len(suc_len)], inverse_data$suc_slacks))
  }

  ## Dual variables for the original equality constraints.
  y_len <- ifelse(length(inverse_data$y_slacks) == 0, 0, sum(sapply(inverse_data$y_slacks, function(val) { val[[2]] })))
  if(y_len > 0) {
      ##y <- rep(0, y_len)
      ## task.getyslice(sol, suc_len, suc_len + y_len, y)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$suc[seq.int(suc_len, length.out = y_len)], inverse_data$y_slacks))
  }

  ## Dual variables for SOC and EXP constraints.
  snx_len <- ifelse(length(inverse_data$snx_slacks) == 0, 0, sum(sapply(inverse_data$snx_slacks, function(val) { val[[2]] })))
  if(snx_len > 0) {
      ##snx <- matrix(0, nrow = snx_len, ncol = 1)
      ##task.getsnxslice(sol, inverse_data$n0, inverse_data$n0 + snx_len, snx)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$snx, inverse_data$snx_slacks))
  }

  ## Dual variables for PSD constraints.
  for(psd_info in inverse_data$psd_dims) {
    id <- as.character(psd_info[[1L]])
    dim <- psd_info[[2L]]
    ##sj <- rep(0, dim*floor((dim + 1)/2))
    ##task.getbars(sol, j, sj)
    dual_vars[[id]] <- vectorized_lower_tri_to_mat(sol$bars[[1L]], dim)
  }

  return(dual_vars)
}

#'
#' Parses MOSEK dual variables into corresponding CVXR constraints and dual values
#'
#' @param dual_var List of the dual variables returned by the MOSEK solution.
#' @param constr_id_to_constr_dim A list that contains the mapping of entry "id"
#' that is the index of the CVXR \linkS4class{Constraint} object to which the
#' next "dim" entries of the dual variable belong.
#' @return A list with the mapping of the CVXR \linkS4class{Constraint} object
#' indices with the corresponding dual values.
MOSEK.parse_dual_vars <- function(dual_var, constr_id_to_constr_dim) {
  dual_vars <- list()
  running_idx <- 1
  for(val in constr_id_to_constr_dim) {
    id <- as.character(val[[1]])
    dim <- val[[2]]
    ## if(dim == 1)
    ##   dual_vars[id] <- dual_vars[running_idx]   # a scalar.
    ## else
    ##   dual_vars[id] <- as.matrix(dual_vars[running_idx:(running_idx + dim)])
    dual_vars[[id]] <- dual_var[seq.int(running_idx, length.out = dim)]
    running_idx <- running_idx + dim
  }
  return(dual_vars)
}

## MOSEK._handle_mosek_params <- function(task, params) {
##     if(is.na(params))
##         return()

##   ##requireNamespace("Rmosek", quietly = TRUE)

##   handle_str_param <- function(param, value) {
##     if(startsWith(param, "MSK_DPAR_"))
##       task.putnadourparam(param, value)
##     else if(startsWith(param, "MSK_IPAR_"))
##       task.putnaintparam(param, value)
##     else if(startsWith(param, "MSK_SPAR_"))
##       task.putnastrparam(param, value)
##     else
##       stop("Invalid MOSEK parameter ", param)
##   }

##   handle_enum_param <- function(param, value) {
##     if(is(param, "dparam"))
##       task.putdouparam(param, value)
##     else if(is(param, "iparam"))
##       task.putintparam(param, value)
##     else if(is(param, "sparam"))
##       task.putstrparam(param, value)
##     else
##       stop("Invalid MOSEK parameter ", param)
##   }

##   for(p in params) {
##     param <- p[[1]]
##     value <- p[[2]]
##     if(is(param, "character"))
##       handle_str_param(param, value)
##     else
##       handle_enum_param(param, value)
##   }
## }

#'
#' Utility methods for special handling of semidefinite constraints.
#'
#' @param matrix The matrix to get the lower triangular matrix for
#' @return The lower triangular part of the matrix, stacked in column-major order
scaled_lower_tri <- function(matrix) {
  # Returns an expression representing the lower triangular entries.
  # Scales the strictly lower triangular entries by sqrt(2), as
  # required by SCS.
  rows <- cols <- nrow(matrix)
  entries <- floor(rows * (cols + 1)/2)

  row_arr <- seq_len(entries)

  col_arr <- matrix(1:(rows*cols), nrow = rows, ncol = cols)
  col_arr <- col_arr[lower.tri(col_arr, diag = TRUE)]

  val_arr <- matrix(0, nrow = rows, ncol = cols)
  val_arr[lower.tri(val_arr, diag = TRUE)] <- sqrt(2)
  diag(val_arr) <- 1
  val_arr <- as.vector(val_arr)
  val_arr <- val_arr[val_arr != 0]

  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(entries, rows*cols)))
  vectorized_matrix <- reshape_expr(matrix, c(rows*cols, 1))
  return(coeff %*% vectorized_matrix)
}

