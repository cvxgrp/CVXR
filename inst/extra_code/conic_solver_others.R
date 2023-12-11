#' An interface for the MOSEK solver.
#'
#' @name MOSEK-class
#' @aliases MOSEK
#' @rdname MOSEK-class
#' @export
setClass("MOSEK", representation(exp_cone_order = "numeric"),   # Order of exponential cone constraints. Internal only!
                  prototype(exp_cone_order = c(2, 1, 0), MIP_CAPABLE = TRUE,
                            SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC", "PSDConstraint", "ExpCone"),
                            MI_SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC")),   # Does not support MISDP.
         contains = "ConicSolver")

#' @rdname MOSEK-class
#' @export
MOSEK <- function() { new("MOSEK") }

#'
#' Turns symmetric 2D array into a lower triangular matrix
#'
#' @param v A list of length (dim * (dim + 1) / 2).
#' @param dim The number of rows (equivalently, columns) in the output array.
#' @return Return the symmetric 2D array defined by taking "v" to specify its
#' lower triangular matrix.
vectorized_lower_tri_to_mat <- function(v, dim) {
  v <- unlist(v)
  rows <- c()
  cols <- c()
  vals <- c()
  running_idx <- 1
  for(j in seq_len(dim)) {
    rows <- c(rows, j + seq_len(dim-j+1) - 1)
    cols <- c(cols, rep(j, dim-j+1))
    vals <- c(vals, v[running_idx:(running_idx + dim - j)])
    running_idx <- running_idx + dim - j + 1
  }
  A <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(dim, dim))
  d <- diag(diag(A))
  A <- A + t(A) - d
  return(A)
}

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

#' @param object,x A \linkS4class{MOSEK} object.
#' @describeIn MOSEK Imports the solver.
#' @importFrom utils packageDescription
setMethod("import_solver", "MOSEK", function(solver) {
    requireNamespace("Rmosek", quietly = TRUE) &&
        (!is.null(utils::packageDescription("Rmosek")$Configured.MSK_VERSION))
    ## TODO: Add exponential cone support.
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
  for(constr in problem@constraints) {
    if(!inherits(constr, supported_constraints(object)))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

#' @param constraints A list of \linkS4class{Constraint} objects for which coefficient
#' andd offset data ("G", "h" respectively) is needed
#' @param exp_cone_order A parameter that is only used when a \linkS4class{Constraint} object
#' describes membership in the exponential cone.
#' @describeIn MOSEK Returns a large matrix "coeff" and a vector of constants "offset" such
#' that every \linkS4class{Constraint} in "constraints" holds at z in R^n iff
#' "coeff" * z <=_K offset", where K is a product of cones supported by MOSEK
#' and CVXR (zero cone, nonnegative orthant, second order cone, exponential cone). The
#' nature of K is inferred later by accessing the data in "lengths" and "ids".
setMethod("block_format", "MOSEK", function(object, problem, constraints, exp_cone_order = NA) {
  if(length(constraints) == 0 || is.null(constraints) || any(is.na(constraints)))
    return(list(NULL, NULL))
  matrices <- list()
  offsets <- c()
  lengths <- c()
  ids <- c()

  for(con in constraints) {
    coeff_offs <- reduction_format_constr(object, problem, con, exp_cone_order)
    coeff <- coeff_offs[[1]]
    offset <- coeff_offs[[2]]
    matrices <- c(matrices, list(coeff))
    offsets <- c(offsets, offset)
    lengths <- c(lengths, prod(dim(as.matrix(offset))))
    ids <- c(ids, id(con))
  }
  coeff <- Matrix(do.call(rbind, matrices), sparse = TRUE)
  return(list(coeff, offsets, lengths, ids))
})

#' @describeIn MOSEK Returns a new problem and data for inverting the new solution.
setMethod("perform", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list(suc_slacks = list(), y_slacks = list(), snx_slacks = list(), psd_dims = list())
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])

  # Get integrality constraint information.
  var <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- sapply(var@boolean_idx, function(t) { as.integer(t[1]) })
  data[[INT_IDX]] <- sapply(var@integer_idx, function(t) { as.integer(t[1]) })
  inv_data$integer_variables <- length(data[[BOOL_IDX]]) + length(data[[INT_IDX]]) > 0

  # Parse the coefficient vector from the objective.
  coeff_offs <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  c <- coeff_offs[[1]]
  constant <- coeff_offs[[2]]
  data[[C_KEY]] <- as.vector(c)
  inv_data$n0 <- length(data[[C_KEY]])
  data[[OBJ_OFFSET]] <- constant[1]
  data[[DIMS]] <- list()
  data[[DIMS]][[SOC_DIM]] <- list()
  data[[DIMS]][[EXP_DIM]] <- list()
  data[[DIMS]][[PSD_DIM]] <- list()
  data[[DIMS]][[LEQ_DIM]] <- 0
  data[[DIMS]][[EQ_DIM]] <- 0
  inv_data[[OBJ_OFFSET]] <- constant[1]
  Gs <- list()
  hs <- list()

  if(length(problem@constraints) == 0) {
    ##data[[G_KEY]] <- Matrix(nrow = 0, ncol = 0, sparse = TRUE)
    ## Ensure G's dimensions match that of c.
    data[[G_KEY]] <- Matrix(nrow = 0, ncol = length(c), sparse = TRUE)
    data[[H_KEY]] <- matrix(nrow = 0, ncol = 0)
    inv_data$is_LP <- TRUE
    return(list(object, data, inv_data))
  }

  # Linear inequalities.
  leq_constr <- problem@constraints[sapply(problem@constraints, inherits, what = "NonPosConstraint" )]
  if(length(leq_constr) > 0) {
    blform <- block_format(object, problem, leq_constr)   # G, h : G*z <= h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$suc_slacks <- c(inv_data$suc_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[[DIMS]][[LEQ_DIM]] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Linear equations.
  eq_constr <- problem@constraints[sapply(problem@constraints, inherits, what = "ZeroConstraint")]
  if(length(eq_constr) > 0) {
    blform <- block_format(object, problem, eq_constr)   # G, h : G*z == h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$y_slacks <- c(inv_data$y_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[[DIMS]][[EQ_DIM]] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Second order cone.
  soc_constr <- problem@constraints[sapply(problem@constraints, inherits, what = "SOC" )]
  data[[DIMS]][[SOC_DIM]] <- list()
  for(ci in soc_constr)
    data[[DIMS]][[SOC_DIM]] <- c(data[[DIMS]][[SOC_DIM]], cone_sizes(ci))
  if(length(soc_constr) > 0) {
    blform <- block_format(object, problem, soc_constr)   # G*z <=_{soc} h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$snx_slacks <- c(inv_data$snx_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Exponential cone.
  exp_constr <- problem@constraints[sapply(problem@constraints, inherits, what = "ExpCone" )]
  if(length(exp_constr) > 0) {
    # G*z <=_{EXP} h.
    blform <- block_format(object, problem, exp_constr, object@exp_cone_order)
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    data[[DIMS]][[EXP_DIM]] <- lengths
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # PSD constraints.
  psd_constr <- problem@constraints[sapply(problem@constraints, inherits, what = "PSDConstraint" )]
  if(length(psd_constr) > 0) {
    data[[DIMS]][[PSD_DIM]] <- list()
    for(c in psd_constr) {
      coeff_offs <- psd_coeff_offset(problem, c)
      G_vec <- coeff_offs[[1]]
      h_vec <- coeff_offs[[2]]
      dim <- coeff_offs[[3]]
      inv_data$psd_dims <- c(inv_data$psd_dims, list(list(id(c), dim)))
      data[[DIMS]][[PSD_DIM]] <- c(data[[DIMS]][[PSD_DIM]], list(dim))
      Gs <- c(Gs, G_vec)
      hs <- c(hs, h_vec)
    }
  }

  if(length(Gs) == 0)
    ## data[[G_KEY]] <- Matrix(nrow = 0, ncol = 0, sparse = TRUE)
    ## G is already sparse
    data[[G_KEY]] <- G
  else
    data[[G_KEY]] <- Matrix(do.call(rbind, Gs), sparse = TRUE)
  if(length(hs) == 0)
    data[[H_KEY]] <- matrix(nrow = 0, ncol = 0)
  else
    data[[H_KEY]] <- Matrix(do.call(cbind, hs), sparse = TRUE)
  inv_data$is_LP <- (length(psd_constr) + length(exp_constr) + length(soc_constr)) == 0
  return(list(object, data, inv_data))
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
setMethod("solve_via_data", "MOSEK", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                              num_iter,
                                              solver_opts, solver_cache) {

  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  ## Check if the CVXR standard form has zero variables. If so,
  ## return a trivial solution. This is necessary because MOSEK
  ## will crash if handed a problem with zero variables.

  c <- data[[C_KEY]]

  if (length(c) == 0) {
      res <- list()
      res[[STATUS]] <- OPTIMAL
      res[[PRIMAL]] <- list()
      res[[VALUE]] <- data[[OFFSET]]
      res[[EQ_DUAL]] <- list()
      res[[INEQ_DUAL]] <- list()
      return(res)
  }
  # The following lines recover problem parameters, and define helper constants.
  #
  #   The problem's objective is "min c.T * z".
  #   The problem's constraint set is "G * z <=_K h."
  #   The rows in (G, h) are formatted in order of
  #       (1) linear inequalities,
  #       (2) linear equations,
  #       (3) soc constraints,
  #       (4) exponential cone constraints,
  #       (5) vectorized linear matrix inequalities.
  #   The parameter "dims" indicates the exact
  #   dimensions of each of these cones.
  #
  #   MOSEK's standard form requires that we replace generalized
  #   inequalities with slack variables and linear equations.
  #   The parameter "n" is the size of the column-vector variable
  #   after adding slacks for SOC and EXP constraints. To be
  #   consistent with MOSEK documentation, subsequent comments
  #   refer to this variable as "x".

  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  dims <- data[[DIMS]]
  dims_SOC_DIM <- dims[[SOC_DIM]]
  dims_EXP_DIM <- dims[[EXP_DIM]]
  dims_PSD_DIM <- dims[[PSD_DIM]]
  data_BOOL_IDX <- data[[BOOL_IDX]]
  data_INT_IDX <- data[[INT_IDX]]
  num_bool <- length(data_BOOL_IDX)
  num_int <- length(data_INT_IDX)
  length_dims_SOC_DIM <- length(dims_SOC_DIM)
  unlist_dims_EXP_DIM <- unlist(dims_EXP_DIM)
  unlist_dims_SOC_DIM <- unlist(dims_SOC_DIM)
  dims_LEQ_DIM <- dims[[LEQ_DIM]]

  n0 <- length(c)
  n <- n0 + sum(unlist_dims_SOC_DIM, na.rm = TRUE) + sum(unlist_dims_EXP_DIM, na.rm = TRUE) # unlisted dims to make sure sum function works and na.rm to handle empty lists
  psd_total_dims <- sum(unlist(dims_PSD_DIM)^2, na.rm = TRUE)
  m <- length(h)



  # Define variables, cone constraints, and integrality constraints.
  #
  #   The variable "x" is a length-n block vector, with
  #       Block 1: "z" from "G * z <=_K h",
  #       Block 2: slacks for SOC constraints, and
  #       Block 3: slacks for EXP cone constraints.
  #
  #   Once we declare x in the MOSEK model, we add the necessary
  #   conic constraints for slack variables (Blocks 2 and 3).
  #   The last step is to add integrality constraints.
  #
  #   Note that the API call for PSD variables contains the word "bar".
  #   MOSEK documentation consistently uses "bar" as a sort of flag,
  #   indicating that a function deals with PSD variables.



  ##env <- Rmosek::Env() remove these as Rmosek doesn't need environments
  ##task <- env.Task(0,0)
  ##instead defines prob
  prob <- list(sense="min")

  ## TODO: Handle logging for verbose.

  ## Parse all user-specified parameters (override default logging
  ## parameters if applicable).
  ## Rmosek expects a list of lists
  ## prob$dparam <- list(...); prob$iparam <- list(...); prob$sparam <- list(...)
  if (!is.null(solver_opts)) {
      prob$dparam  <-  solver_opts$dparam
      prob$iparam  <-  solver_opts$iparam
      prob$sparam  <-  solver_opts$sparam
  }

  if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol), is.null(num_iter)))) {
    warning("Ignoring inapplicable parameter feastol/reltol/abstol/num_iter for MOSEK.")
  }


  #task.appendvars(n), no task for Rmosek, but declares the number of variables in the model. Need to expand prob$c as well to match this dimension
  #task.putvarboundlist(1:n, rep(mosek.boundkey.fr, n), matrix(0, nrow = n, ncol = 1), matrix(0, nrow = n, ncol = 1))
  #kind of confused why x's are all 0's, but that's what's in python code
  #prob$bx <- rbind( blx = rep(0,n),
  #                 bux = rep(0,n))
  prob$bx <- rbind(blx = rep(-Inf, n),
                   bux = rep(Inf, n))

  #Initialize the cone. Not 100% sure about this bit
  NUMCONES <- length_dims_SOC_DIM + floor(sum(unlist_dims_EXP_DIM, na.rm = TRUE)/3)
  prob$cones <- matrix(list(), nrow = 2, ncol = NUMCONES)

  if(psd_total_dims > 0)
    prob$bardim <- unlist(dims_PSD_DIM)
  running_idx <- n0
  for(i in seq_along(unlist_dims_SOC_DIM)) {
    prob$cones[,i] <- list("QUAD", as.numeric((running_idx + 1):(running_idx + unlist_dims_SOC_DIM[[i]]))) # latter term is size_cone
    running_idx <- running_idx + unlist_dims_SOC_DIM[[i]]
  }
  if(floor(sum(unlist_dims_EXP_DIM, na.rm = TRUE)/3) != 0){ # check this, feels sketchy
    for(k in 1:floor(sum(unlist_dims_EXP_DIM, na.rm = TRUE)/3) ) {
      prob$cones[,(length_dims_SOC_DIM+k)] <- list("PEXP", as.numeric((running_idx+1):(running_idx + 3)) )
      running_idx <- running_idx + 3
    }
  }
  if(num_bool + num_int > 0) {
    if(num_bool > 0) {
      unlist_data_BOOL_IDX <- unlist(data_BOOL_IDX)
      prob$intsub <- unlist_data_BOOL_IDX
      #since the variable constraints are already declared, we are resetting them so they can only be 0 or 1
      prob$bx[, unlist_data_BOOL_IDX] <- rbind( rep(0, length(unlist_data_BOOL_IDX)), rep(1, length(unlist_data_BOOL_IDX)) )

    }
    if(num_int > 0)
      prob$intsub <- unlist(data_INT_IDX)
  }

  # Define linear inequality and equality constraints.
  #
  #   Mosek will see a total of m linear expressions, which must
  #   define linear inequalities and equalities. The variable x
  #   contributes to these linear expressions by standard
  #   matrix-vector multiplication; the matrix in question is
  #   referred to as "A" in the mosek documentation. The PSD
  #   variables have a different means of contributing to the
  #   linear expressions. Specifically, a PSD variable Xj contributes
  #   "+tr( \bar{A}_{ij} * Xj )" to the i-th linear expression,
  #   where \bar{A}_{ij} is specified by a call to putbaraij.
  #
  #   The following code has three phases.
  #       (1) Build the matrix A.
  #       (2) Specify the \bar{A}_{ij} for PSD variables.
  #       (3) Specify the RHS of the m linear (in)equalities.
  #
  #   Remark : The parameter G gives every row in the first
  #   n0 columns of A. The remaining columns of A are for SOC
  #   and EXP slack variables. We can actually account for all
  #   of these slack variables at once by specifying a giant
  #   identity matrix in the appropriate position in A.

  # task.appendcons(m) is equivalent to prob$bc


  ##G should already be sparse but Matrix 1.3.x causes problems.
  ## if (!inherits(G, "dgCMatrix")) G  <- as(as(G, "CsparseMatrix"), "dgCMatrix")
  ## Matrix 1.5 change
  if (!inherits(G, "dgCMatrix")) G  <- as(as(G, "CsparseMatrix"), "generalMatrix")

  G_sum <- summary(G)
  nrow_G_sparse <- nrow(G)
  ncol_G_sparse <- ncol(G)

  row <- G_sum$i
  col <- G_sum$j
  vals <- G_sum$x

  total_soc_exp_slacks <- sum(unlist_dims_SOC_DIM, na.rm = TRUE) + sum(unlist_dims_EXP_DIM, na.rm = TRUE)

  # initializing A matrix
  if(ncol_G_sparse == 0 || (ncol_G_sparse + total_soc_exp_slacks) == 0)
    ## prob$A <- sparseMatrix(i = c(), j = c(), dims = c(0, 0))
    ## G is already sparse
    prob$A  <- G
  else {
    # this is a bit hacky, probably should fix later. Filling out part of the A matrix from G
    # Equivalent to task.putaijlist(as.list(row), as.list(col), as.list(vals))
    if(total_soc_exp_slacks > 0) {
      i <- unlist(dims[[LEQ_DIM]]) + unlist(dims[[EQ_DIM]])   # Constraint index in (1, ..., m)
      j <- length(c)   # Index of the first slack variable in the block vector "x".
      rows <- (i:(i + total_soc_exp_slacks-1))+1
      cols <- (j:(j + total_soc_exp_slacks-1))+1

      row <- c(row, rows)
      col <- c(col, cols)
      vals <- c(vals, rep(1, length(rows)))
    }
    prob$A <- sparseMatrix(i = row, j = col, x = vals, dims = c(nrow_G_sparse, ncol_G_sparse + total_soc_exp_slacks))
  }

  # Constraint index: start of LMIs.
  i <- dims_LEQ_DIM + dims[[EQ_DIM]] + total_soc_exp_slacks + 1
  dim_exist_PSD <- length(dims_PSD_DIM) #indicates whether or not we have any LMIs

  if(dim_exist_PSD > 0){
    #A bit hacky here too, specifying the lower triangular part of symmetric coefficient matrix barA
    barAi <- c() #Specifies row index of block matrix
    barAj <- c() #Specifies column index of block matrix
    barAk <- c() #Specifies row index within the block matrix specified above
    barAl <- c() #Specifies column index within the block matrix specified above
    barAv <- c() #Values for all the matrices

    for(j in 1:length(dims_PSD_DIM)) {   #For each PSD matrix
      for(row_idx in 1:dims_PSD_DIM[[j]]) {
        for(col_idx in 1:dims_PSD_DIM[[j]]) {
          val <- ifelse(row_idx == col_idx, 1, 0.5)
          row <- max(row_idx, col_idx)
          col <- min(row_idx, col_idx)
          #mat <- task.appendsparsesymmat(dim, list(row), list(col), list(val))
          #task.putbaraij(i, j, list(mat), list(1.0))

          barAi <- c(barAi, i)
          barAj <- c(barAj, j) #NEED TO CHECK. Multiple PSD_DIM example?
          barAk <- c(barAk, row)
          barAl <- c(barAl, col)
          barAv <- c(barAv, val)

          i <- i + 1 #for each symmetric matrix
        }
      }
    }

    #Attaching. Does mosek automatically check the symmetric matrix dimensions?

    prob$barA$i <- barAi
    prob$barA$j <- barAj
    prob$barA$k <- barAk
    prob$barA$l <- barAl
    prob$barA$v <- barAv
  }

  num_eq <- length(h) - dims_LEQ_DIM

  #CVXPY has the first dims[[LEQ_DIM]] variables as upper bounded
  #type_constraint <- rep(mosek.boundkey.up, dims[[LEQ_DIM]]) + rep(mosek.boundkey.fx, num_eq)
  #task.putconboundlist(1:m, type_constraint, h, h), equivalent to prob$bc
  hl_holder <- as.numeric(h)
  hu_holder <- as.numeric(h)

  #upper constraints for the LEQ_DIM number of variables, so set lower bound to -Inf
  hl_holder[seq_len(dims_LEQ_DIM)] <- rep(-Inf, dims_LEQ_DIM)

  prob$bc <- rbind(blc = hl_holder,
                   buc = hu_holder)

  # Define the objective and optimize the MOSEK task.
  #initialize coefficients of objective with the same number of variables declared (dim of x)
  c_holder <- rep(0, n)
  c_holder[1:length(c)] <- c

  prob$c <- c_holder

  if(is.logical(verbose) && verbose ){
    verbose <- 10
  } else if(!verbose){
    verbose <- 0
  } else if(!is.null(solver_opts$verbose)){
    verbose <- solver_opts$verbose
  }

  if(is.null(solver_opts$soldetail)){
    solver_opts$soldetail <- 3
  } else {
    warning("Solver might not output correct answer depending on the input of the soldetail variable. Default is 3")
  }

  if(is.null(solver_opts$getinfo)){
    solver_opts$getinfo <- TRUE
  } else {
    warning("Solver might not output correct answer depending on the input of the getinfo variable. Default is TRUE")
  }

  r <- Rmosek::mosek(prob, list(verbose = verbose, usesol = solver_opts$usesol,
                                useparam = solver_opts$useparam, soldetail = solver_opts$soldetail,
                                getinfo = solver_opts$getinfo, writebefore = solver_opts$writebefore,
                                writeafter = solver_opts$writeafter))

  return(r)
})

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

#' Utility method for formatting a ConeDims instance into a dictionary
#' that can be supplied to SCS.
#' @param cone_dims A \linkS4class{ConeDims} instance.
#' @return The dimensions of the cones.
SCS.dims_to_solver_dict <- function(cone_dims) {
  cones <- list(z = as.integer(cone_dims@zero),
                l = as.integer(cone_dims@nonpos),
                q = sapply(cone_dims@soc, as.integer),
                ep = as.integer(cone_dims@exp),
                s = sapply(cone_dims@psd, as.integer))
  return(cones)
}

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

# SuperSCS <- setClass("SuperSCS", contains = "SCS")
# SuperSCS.default_settings <- function(object) {
#   list(use_indirect = FALSE, eps = 1e-8, max_iters = 10000)
# }
#
# #' @param solver,object,x A \linkS4class{SuperSCS} object.
# #' @describeIn SuperSCS Returns the name of the SuperSCS solver
# setMethod("name", "SuperSCS", function(x) { SUPER_SCS_NAME })
#
# #' @describeIn SuperSCS imports the SuperSCS solver.
# setMethod("import_solver", "SuperSCS", function(solver) {
#   stop("Unimplemented: SuperSCS is currently unavailable in R.")
# })
#
# #' @param data Data generated via an apply call.
# #' @param warm_start An option for warm start.
# #' @param verbose A boolean of whether to enable solver verbosity.
# #' @param solver_opts A list of Solver specific options
# #' @param solver_cache Cache for the solver.
# #' @describeIn SuperSCS Solve a problem represented by data returned from apply.
# setMethod("solve_via_data", "SuperSCS", function(object, data, warm_start, verbose, solver_opts, solver_cache) {
#   if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
#   args <- list(A = data[[A_KEY]], b = data[[B_KEY]], c = data[[C_KEY]])
#   if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
#     args$x <- solver_cache[[name(object)]]$x
#     args$y <- solver_cache[[name(object)]]$y
#     args$s <- solver_cache[[name(object)]]$s
#   }
#   cones <- SCS.dims_to_solver_dict(data[[ConicSolver()@dims]])
#
#   # Settings.
#   user_opts <- names(solver_opts)
#   for(k in names(SuperSCS.default_settings)) {
#     if(!k %in% user_opts)
#       solver_opts[[k]] <- SuperSCS.default_settings[[k]]
#   }
#   results <- SuperSCS::solve(args, cones, verbose = verbose, solver_opts)
#   if(!is.null(solver_cache) && length(solver_cache) > 0)
#     solver_cache[[name(object)]] <- results
#   return(results)
# })

# XPRESS <- setClass("XPRESS", prototype(MIP_CAPABLE = TRUE, SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC"),
#                                        MI_SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC")), contains = "SCS")
#
# # Map of XPRESS status to CVXR status.
# setMethod("status_map", "XPRESS", function(solver, status) {
#   if(status == 2)
#     return(OPTIMAL)
#   else if(status == 3)
#     return(INFEASIBLE)
#   else if(status == 5)
#     return(UNBOUNDED)
#   else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
#     return(SOLVER_ERROR)
#   else if(status == 9)   # TODO: Could be anything. Means time expired.
#     return(OPTIMAL_INACCURATE)
#   else
#     stop("XPRESS status unrecognized: ", status)
# })
#
# setMethod("name", "XPRESS", function(x) { XPRESS_NAME })
# setMethod("import_solver", "XPRESS", function(solver) {
#   stop("Unimplemented: XPRESS solver unavailable in R.")
# })
#
# setMethod("accepts", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
#   # TODO: Check if the matrix is stuffed.
#   if(!is_affine(problem@objective@args[[1]]))
#     return(FALSE)
#   for(constr in problem@constraints) {
#     if(!inherits(constr, supported_constraints(object)))
#       return(FALSE)
#     for(arg in constr@args) {
#       if(!is_affine(arg))
#         return(FALSE)
#     }
#   }
#   return(TRUE)
# })
#
# setMethod("perform", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
#   tmp <- callNextMethod(object, problem)
#   data <- tmp[[1]]
#   inv_data <- tmp[[2]]
#   variables <- variables(problem)[[1]]
#   data[[BOOL_IDX]] <- lapply(variables@boolean_idx, function(t) { t[1] })
#   data[[INT_IDX]] <- lapply(variables@integer_idx, function(t) { t[1] })
#   inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0
#   return(list(object, data, inv_data))
# })
#
# setMethod("invert", signature(object = "XPRESS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
#   status <- solution[[STATUS]]
#
#   if(status %in% SOLUTION_PRESENT) {
#     opt_val <- solution[[VALUE]]
#     primal_vars <- list()
#     primal_vars[[inverse_data[[object@var_id]]]] <- solution$primal
#     if(!inverse_data@is_mip)
#       dual_vars <- get_dual_values(solution[[EQ_DUAL]], extract_dual_value, inverse_data[[EQ_CONSTR]])
#   } else {
#     primal_vars <- list()
#     primal_vars[[inverse_data[[object@var_id]]]] <- NA_real_
#     if(!inverse_data@is_mip) {
#       dual_var_ids <- sapply(inverse_data[[EQ_CONSTR]], function(constr) { constr@id })
#       dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
#       names(dual_vars) <- dual_var_ids
#     }
#
#     if(status == INFEASIBLE)
#       opt_val <- Inf
#     else if(status == UNBOUNDED)
#       opt_val <- -Inf
#     else
#       opt_val <- NA
#   }
#
#   other <- list()
#   other[[XPRESS_IIS]] <- solution[[XPRESS_IIS]]
#   other[[XPRESS_TROW]] <- solution[[XPRESS_TROW]]
#   return(Solution(status, opt_val, primal_vars, dual_vars, other))
# })
#
# setMethod("solve_via_data", "XPRESS", function(object, data, warm_start, verbose, solver_opts, solver_cache) {
#   if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
#   solver <- XPRESS_OLD()
#   solver_opts[[BOOL_IDX]] <- data[[BOOL_IDX]]
#   solver_opts[[INT_IDX]] <- data[[INT_IDX]]
#   prob_data <- list()
#   prob_data[[name(object)]] <- ProblemData()
#   solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
# })
