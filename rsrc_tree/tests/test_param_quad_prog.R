context("test_param_quad_prog")
TOL <- 1e-2

solvers <- installed_solvers()
solvers <- solvers[solvers %in% CVXR:::QP_SOLVERS]

test_that("test param data", {
  for(solver in solvers) {
    set.seed(0)
    m <- 30
    n <- 20
    A <- matrix(rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(rnorm(m), nrow = m, ncol = 1)
    x <- Variable(n)
    gamma <- Parameter(nonneg = TRUE)
    gamma_val <- 0.5
    gamma_val_new <- 0.1
    objective <- Minimize(gamma*sum_squares(A %*% x - b) + norm1(x))
    constraints <- list(1 <= x, x <= 2)
    
    # Solve from scratch (directly new parameter).
    prob <- Problem(objective, constraints)
    expect_true(is_dpp(prob))
    value(gamma) <- gamma_val_new
    data_scratch <- get_problem_data(prob, solver)[[1]]
    result <- solve(prob, solver = solver)
    x_scratch <- value(x)
    
    # Canonicalize problem with parameter values (solve once).
    prob <- Problem(objective, constraints)
    value(gamma) <- gamma_val
    data_param <- get_problem_data(prob, solver)[[1]]
    result <- solve(prob, solver = solver)
    
    # Get data with new parameter.
    value(gamma) <- gamma_val_new
    data_param_new <- get_problem_data(prob, solver)[[1]]
    result <- solve(prob, solver = solver)
    x_gamma_new <- value(x)
    
    # Check if data match.
    expect_true(is.allclose(as.matrix(data_param_new$P), as.matrix(data_scratch$P)))
    
    # Check if solutions match.
    expect_true(is.allclose(x_gamma_new, x_scratch, rtol = 1e-2, atol = 1e-2))
  }
})

test_that("test qp problem", {
  for(solver in solvers) {
    m <- 30
    n <- 20
    A <- matrix(rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(rnorm(m), nrow = m, ncol = 1)
    x <- Variable(n)
    gamma <- Parameter(nonneg = TRUE)
    value(gamma) <- 0.5
    objective <- Minimize(sum_squares(A %*% x - b) + gamma*norm1(x))
    constraints <- list(0 <= x, x <= 1)
    
    # Solve from scratch.
    problem <- Problem(objective, constraints)
    result <- solve(problem, solver = solver)
    x_full <- value(x)
    
    # Restore cached values.
    solving_chain <- problem@.cache@solving_chain
    solver <- problem@.cache@solving_chain@solver
    inverse_data <- problem@.cache@inverse_data
    param_prog <- problem@.cache@param_prog
    
    # Solve parameteric
    tmp <- perform(solving_chain@solver, param_prog)
    solving_chain@solver <- tmp[[1]]
    data <- tmp[[2]]
    solver_inverse_data <- tmp[[3]]
    inverse_data <- c(inverse_data, list(solver_inverse_data))
    raw_solution <- solve_via_data(solver, data, warm_start = FALSE, verbose = FALSE, solver_opts = list())
    result <- unpack_results(problem, raw_solution, solving_chain, inverse_data)
    x_param <- result$getValue(x)
    
    expect_true(is.allclose(x_param, x_full, rtol = 1e-2, atol = 1e-2))
  }
  
  # TODO: Add derivatives and adjoint tests.
})
