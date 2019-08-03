context("test-g01-solvers")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test that all the ECOS solver options work", {
  # Test ecos
  # feastol, abstol, reltol, feastol_inacc,
  # abstol_inacc, and reltol_inacc for tolerance values
  # max_iters for the maximum number of iterations,
  EPS <- 1e-4
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
  for(i in 1:2) {
    result <- solve(prob, solver = "ECOS", feastol = EPS, abstol = EPS, reltol = EPS, 
                    feastol_inacc = EPS, abstol_inacc = EPS, reltol_inacc = EPS, 
                    max_iters = 20, verbose = TRUE, warm_start = TRUE)
  }
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
})

test_that("test that all the ECOS BB solver options work", {
  # 'mi_maxiter'
  # maximum number of branch and bound iterations (default: 1000)
  # 'mi_abs_eps'
  # absolute tolerance between upper and lower bounds (default: 1e-6)
  # 'mi_rel_eps'
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == Variable(2, boolean = TRUE)))
  for(i in 1:2) {
    result <- solve(prob, solver = "ECOS_BB", mi_max_iters = 100, mi_abs_eps = 1e-6,
                    mi_rel_eps = 1e-5, verbose = TRUE, warm_start = TRUE)
  }
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
})

test_that("test that all the SCS solver options work", {
  # Test SCS
  # MAX_ITERS, EPS, ALPHA, UNDET_TOL, VERBOSE, and NORMALIZE.
  # If opts is missing, then the algorithm uses default settings.
  # USE_INDIRECT = True
  EPS <- 1e-4
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
  for(i in 1:2) {
    result <- solve(prob, solver = "SCS", max_iters = 50, eps = EPS, alpha = EPS,
                    verbose = TRUE, normalize = TRUE, use_indirect = FALSE)
  }
  expect_equal(result$value, 1.0, tolerance = 1e-2)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = 1e-2)
})

test_that("test a basic LP with GUROBI", {
  if("GUROBI" %in% installed_solvers()) {
    prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
    
    # Gurobi's default lower bound for a decision variable is zero
    # This quick test ensures that the cvxpy interface for GUROBI does *not* have that bound
    objective <- Minimize(x[1])
    constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
    
    # Boolean and integer version.
    bool_var <- Variable(boolean = TRUE)
    int_var <- Variable(integer = TRUE)
    prob <- Problem(Minimize(p_norm(x,1)), list(x == bool_var, bool_var == 0))
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= int_var, x[1] + 2*x[2] <= 3*bool_var, x[1] >= 0, x[2] >= 0, int_var == 3*bool_var, int_var == 3)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(int_var), 3, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "GUROBI"), "The solver GUROBI is not installed.")
  }
})

test_that("test a basic SOCP with GUROBI", {
  if("GUROBI" %in% installed_solvers()) {
    prob <- Problem(Minimize(p_norm(x,2) + 1.0), list(x == 0))
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
    
    # GUROBI's default lower bound for a decision variable is zero
    # This quick test ensures that the CVXR interface for GUROBI does *not* have that bound
    objective <- Minimize(x[1])
    constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
    
    # Boolean and integer version.
    bool_var <- Variable(boolean = TRUE)
    int_var <- Variable(integer = TRUE)
    prob <- Problem(Minimize(p_norm(x,2)), list(x == bool_var, bool_var == 0))
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, 0, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= int_var, (x[1] + 2*x[2])^2 <= 9*bool_var, x[1] >= 0, x[2] >= 0, int_var == 3*bool_var, int_var == 3)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(int_var), 3, tolerance = TOL)
    expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "GUROBI"), "The solver GUROBI is not installed.")
  }
})

test_that("Make sure GUROBI's dual result matches other solvers", {
  if("GUROBI" %in% installed_solvers()) {
    constraints <- list(x == 0)
    prob <- Problem(Minimize(p_norm(x,1)))
    result <- solve(prob, solver = "GUROBI")
    duals_gurobi <- lapply(constraints, function(c) { result$getDualValue(c) })
    result <- solve(prob, solver = "ECOS")
    duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
    for(i in seq_along(constraints))
      expect_equal(duals_gurobi[[i]], duals_ecos[[i]], tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "GUROBI")
    duals_gurobi <- lapply(constraints, function(c) { result$getDualValue(c) })
    result <- solve(prob, solver = "ECOS")
    duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
    for(i in seq_along(constraints))
      expect_equal(duals_gurobi[[i]], duals_ecos[[i]], tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "GUROBI"), "The solver GUROBI is not installed.")
  }
})

test_that("test a basic LP with MOSEK", {
  if("MOSEK" %in% installed_solvers()) {
    prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
    
    # Gurobi's default lower bound for a decision variable is zero
    # This quick test ensures that the cvxpy interface for GUROBI does *not* have that bound
    objective <- Minimize(x[1])
    constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "MOSEK"), "The solver MOSEK is not installed.")
  }
})

test_that("test a basic SOCP with MOSEK", {
  if("MOSEK" %in% installed_solvers()) {
    prob <- Problem(Minimize(p_norm(x,2) + 1.0), list(x == 0))
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
    
    objective <- Minimize(x[1])
    constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "MOSEK"), "The solver MOSEK is not installed.")
  }
})

test_that("Make sure MOSEK's dual result matches other solvers", {
  if("MOSEK" %in% installed_solvers()) {
    constraints <- list(x == 0)
    prob <- Problem(Minimize(p_norm(x,1)))
    result <- solve(prob, solver = "MOSEK")
    duals_mosek <- lapply(constraints, function(c) { result$getDualValue(c) })
    result <- solve(prob, solver = "ECOS")
    duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
    for(i in seq_along(constraints))
      expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = TOL)
    
    # Example from
    # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
    objective <- Minimize(-4*x[1] - 5*x[2])
    constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
    prob <- Problem(objective, constraints)
    result <- solve(prob, solver = "MOSEK")
    duals_mosek <- lapply(constraints, function(c) { result$getDualValue(c) })
    result <- solve(prob, solver = "ECOS")
    duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
    for(i in seq_along(constraints))
      expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = TOL)
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "MOSEK"), "The solver MOSEK is not installed.")
  }
})

test_that("Test a basic SDP with MOSEK", {
  # TODO: Should work with PSD (>>, <<).
  if("MOSEK" %in% installed_solvers()) {
    # Test optimality gap for equilibration.
    n <- 3
    Art <- matrix(rnorm(n*n), nrow = n, ncol = n)
    
    t <- Variable()
    d <- Variable(n)
    D <- diag(d)
    constr <- list(Art %*% D %*% t(Art) - diag(n) == Variable(n, n, PSD = TRUE),
                   Variable(n, n, PSD = TRUE) == t*diag(n) - Art %*% D %*% t(Art), d >= 0)
    prob <- Problem(Minimize(t), constr)
    result <- solve(prob, solver = "MOSEK")
    expect_equal(result$status, "optimal")
  } else {
    prob <- Problem(Minimize(p_norm(x,1)), list(x == 0))
    expect_error(result <- solve(prob, solver = "MOSEK"), "The solver MOSEK is not installed.")
  }
})

