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

SOLVER_MAP_CONIC <- CVXR:::SOLVER_MAP_CONIC
SOLVER_MAP_QP <- CVXR:::SOLVER_MAP_QP
INSTALLED_SOLVERS <- installed_solvers()
CVXOPT_INSTALLED  <- "CVXOPT" %in% INSTALLED_SOLVERS
GLPK_INSTALLED  <- "GLPK"  %in% INSTALLED_SOLVERS
GLPK_MI_INSTALLED  <- "GLPK_MI" %in% INSTALLED_SOLVERS
CPLEX_INSTALLED  <- "CPLEX" %in% INSTALLED_SOLVERS
GUROBI_INSTALLED  <- "GUROBI" %in% INSTALLED_SOLVERS
MOSEK_INSTALLED  <- "MOSEK" %in% INSTALLED_SOLVERS
XPRESS_INSTALLED  <- "XPRESS" %in% INSTALLED_SOLVERS
NAG_INSTALLED  <- "NAG" %in% INSTALLED_SOLVERS

## For CRAN drop CPLEX
##SOLVER_MAP_CONIC$CPLEX  <- NULL
##SOLVER_MAP_QP$CPLEX  <- NULL

test_that("Test that all the ECOS solver options work", {
  skip_on_cran()
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

test_that("Test that all the ECOS BB solver options work", {
  skip_on_cran()
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

test_that("Test that all the SCS solver options work", {
  skip_on_cran()
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

test_that("Test that all the CVXOPT solver options work", {
    skip_on_cran()
    skip_if_not(CVXOPT_INSTALLED)
    prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
    EPS <- 1e-7
    for(i in 1:2) {
      result <- solve(prob, solver = "CVXOPT", feastol = EPS, abstol = EPS, reltol = EPS,
                      max_iters = 20, verbose = TRUE, kktsolver = "chol", refinement = 2, warm_start = TRUE)
    }
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(0, 0)), tolerance = TOL)
})

test_that("Test a basic LP with GLPK", {
  skip_on_cran()
  # Either the problem is solved or GLPK is not installed.
  skip_if_not(GLPK_INSTALLED)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  result <- solve(prob, solver = "GLPK")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(0, 0)), tolerance = TOL)

  ## Example from
  ## http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3,
                      x[1] + 2*x[2] <= 3,
                      x[1] >= 0,
                      x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GLPK")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)
})

test_that("Test a basic MILP with GLPK", {
  skip_on_cran()
  # Either the problem is solved or GLPK is not installed.
  skip_if_not(GLPK_MI_INSTALLED)
  bool_var <- Variable(boolean = TRUE)
  int_var <- Variable(integer = TRUE)
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == bool_var, bool_var == 0))
  result <- solve(prob, solver = "GLPK_MI", verbose = TRUE)
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(0, 0)), tolerance = TOL)

  ## Example from
  ## http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= int_var,
                      x[1] + 2*x[2] <= 3*bool_var,
                      x[1] >= 0,
                      x[2] >= 0,
                      int_var == 3*bool_var,
                      int_var == 3)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GLPK_MI", verbose = TRUE)
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(int_var), 3, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)
})

test_that("Test a basic LP with CPLEX", {
  skip_on_cran()
  skip_if_not(CPLEX_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  ## Example from
  ## http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  ## CPLEX's default lower bound for a decision variable is zero
  ## This quick test ensures that the cvxpy interface for CPLEX does *not* have that bound
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)

  ## Boolean and integer version.
  bool_var <- Variable(boolean = TRUE)
  int_var <- Variable(integer = TRUE)
  prob <- Problem(Minimize(p_norm(x,1)), list(x == bool_var, bool_var == 0))
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  ## Example from
  ## http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= int_var, x[1] + 2*x[2] <= 3*bool_var, x[1] >= 0, x[2] >= 0, int_var == 3*bool_var, int_var == 3)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(int_var), 3, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = TOL)
})

test_that("Test a basic SOCP with CPLEX", {
  skip_on_cran()
  skip_if_not(CPLEX_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,2) + 1.0), list(x == 0))
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  # CPLEX's default lower bound for a decision variable is zero
  # This quick test ensures that the CVXR interface for CPLEX does *not* have that bound
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)

  # Boolean and integer version.
  bool_var <- Variable(boolean = TRUE)
  int_var <- Variable(integer = TRUE)
  prob <- Problem(Minimize(p_norm(x,2)), list(x == bool_var, bool_var == 0))
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= int_var, (x[1] + 2*x[2])^2 <= 9*bool_var, x[1] >= 0, x[2] >= 0, int_var == 3*bool_var, int_var == 3)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(int_var), 3, tolerance = TOL)
  expect_equal(result$getValue(bool_var), 1, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = TOL)
})

test_that("Make sure CPLEX's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(CPLEX_INSTALLED)
  constraints <- list(x == 0)
  prob <- Problem(Minimize(p_norm(x,1)))
  result <- solve(prob, solver = "CPLEX")
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
  result <- solve(prob, solver = "CPLEX")
  duals_gurobi <- lapply(constraints, function(c) { result$getDualValue(c) })
  result <- solve(prob, solver = "ECOS")
  duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
  for(i in seq_along(constraints))
      expect_equal(duals_gurobi[[i]], duals_ecos[[i]], tolerance = TOL)
})

test_that("Test CPLEX warm start", {
  ## Make sure that warm starting CPLEX behaves as expected.
  ## Note: This only checks output, not whether or not CPLEX is warm starting internally.
  skip_on_cran()
  skip_if_not(CPLEX_INSTALLED)
  A <- Parameter(2, 2)
  b <- Parameter(2)
  h <- Parameter(2)
  c <- Parameter(2)

  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(2,2)
  value(c) <- c(1,1)

  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
  expect_equal(result$value, 3)
  expect_equal(result$getValue(x), matrix(c(1, 2)), tolerance = TOL)

  # Change A and b from the original values.
  value(A) <- rbind(c(0,0), c(0,1))   # <----- Changed.
  value(b) <- c(0,1)   # <----- Changed.
  value(h) <- c(2,2)
  value(c) <- c(1,1)

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
  expect_equal(result$value, 3)
  expect_equal(result$getValue(x), matrix(c(2, 1)), tolerance = TOL)

  # Change h from the original values.
  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(1,1)   # <----- Changed.
  value(c) <- c(1,1)

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
  expect_equal(result$value, 2)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)

  # Change c from the original values.
  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(2,2)
  value(c) <- c(2,1)   # <----- Changed.

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "CPLEX", warm_start = TRUE)
  expect_equal(result$value, 4)
  expect_equal(result$getValue(x), matrix(c(1, 2)), tolerance = TOL)
})

## We need a unified interface for solver parameters before doing these checks
## So commenting them out for now

## test_that("Test CPLEX parameters", {
## skip_if_not(CPLEX_INSTALLED)
##   n <- 10
##   m <- 4
##   A <- matrix(rnorm(m*n), nrow = m, ncol = n)
##   x <- matrix(rnorm(n), nrow = n, ncol = 1)
##   y <- A %*% x

##   # Solve a simple basis pursuit problem for testing purposes.
##   z <- Variable(n)
##   objective <- Minimize(norm1(z))
##   constraints <- list(A %*% z == y)
##   problem <- Problem(objective, constraints)
##   ## Until we do a careful check of solver arguments, this check should be commented out
##   ##expect_error(result <- solve(problem, solver = "CPLEX", bogus = "foo"))
##   ##expect_error(result <- solve(problem, solver = "CPLEX", invalid_kwarg = NA))
##   # solve(problem, solver = "CPLEX", advance = 0, simplex.limits.iterations = 1000, timelimit = 1000.0, workdir = "mydir")
## })

test_that("Make sure CVXOPT's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(CVXOPT_INSTALLED)
  constraints <- list(x == 0)
  prob <- Problem(Minimize(p_norm(x,1)))
  result <- solve(prob, solver = "CVXOPT")
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
  result <- solve(prob, solver = "CVXOPT")
  duals_gurobi <- lapply(constraints, function(c) { result$getDualValue(c) })
  result <- solve(prob, solver = "ECOS")
  duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
  for(i in seq_along(constraints))
    expect_equal(duals_gurobi[[i]], duals_ecos[[i]], tolerance = TOL)
})

test_that("Test a basic LP with GUROBI", {
  skip_on_cran()
  skip_if_not(GUROBI_INSTALLED)
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

  # GUROBI's default lower bound for a decision variable is zero
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
})

test_that("Test a basic SOCP with GUROBI", {
  skip_on_cran()
  skip_if_not(GUROBI_INSTALLED)
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
})

test_that("Make sure GUROBI's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(GUROBI_INSTALLED)
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
})

test_that("Test a basic LP with MOSEK", {
  skip_on_cran()
  skip_if_not(MOSEK_INSTALLED)
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

  # MOSEK's default lower bound for a decision variable is zero
  # This quick test ensures that the cvxpy interface for MOSEK does *not* have that bound
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
})

test_that("Test a basic SOCP with MOSEK", {
  skip_on_cran()
  skip_if_not(MOSEK_INSTALLED)
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
})

test_that("Make sure MOSEK's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(MOSEK_INSTALLED)
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
})

test_that("Test a basic SDP with MOSEK", {
  skip_on_cran()
  # TODO: Should work with PSD (>>, <<).
  skip_if_not(MOSEK_INSTALLED)
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
})

## We ignore these test of Mosek parameters since the R mosek interface does not check them.

## test_that("Test MOSEK parameters", {
##   skip_if_not(MOSEK_INSTALLED)
##   n <- 1000
##   m <- 400
##   A <- matrix(rnorm(m*n), nrow = m, ncol = n)
##   x <- matrix(rnorm(n), nrow = n, ncol = 1)
##   y <- A %*% x

##   # Solve a simple basis pursuit problem for testing purposes.
##   z <- Variable(n)
##   objective <- Minimize(norm1(z))
##   constraints <- list(A %*% z == y)
##   problem <- Problem(objective, constraints)
##   ## These tests are currently not useful because the Rmosek interface doesn't check them.
##   ## We have to manually check later
##   ##expect_error(result <- solve(problem, solver = "MOSEK", list(dparam = list(BASIS_TOL_X = "1e-8")))
##   ## expect_error(result <- solve(problem, solver = "MOSEK", list(dparam = list(invalid_kwarg = NA))))
##   solve(problem, solver = "MOSEK",
##         list(dparam = list(BASIS_TOL_X = 1e-1), iparam = list(INTPNT_MAX_ITERATIONS = 20)))
## })

test_that("Test GUROBI warm start", {
  skip_on_cran()
  # Make sure that warm starting GUROBI behaves as expected.
  # Note: This only checks output, not whether or not GUROBI is warm starting internally.
  skip_if_not(GUROBI_INSTALLED)
  A <- Parameter(2, 2)
  b <- Parameter(2)
  h <- Parameter(2)
  c <- Parameter(2)

  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(2,2)
  value(c) <- c(1,1)

  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(result$value, 3)
  expect_equal(result$getValue(x), matrix(c(1, 2)), tolerance = TOL)

  # Change A and b from the original values.
  value(A) <- rbind(c(0,0), c(0,1))   # <----- Changed.
  value(b) <- c(0,1)   # <----- Changed.
  value(h) <- c(2,2)
  value(c) <- c(1,1)

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(result$value, 3)
  expect_equal(result$getValue(x), matrix(c(2, 1)), tolerance = TOL)

  # Change h from the original values.
  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(1,1)   # <----- Changed.
  value(c) <- c(1,1)

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(result$value, 2)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)

  # Change c from the original values.
  value(A) <- rbind(c(1,0), c(0,0))
  value(b) <- c(1,0)
  value(h) <- c(2,2)
  value(c) <- c(2,1)   # <----- Changed.

  # Without setting update_eq_constrs = FALSE, the results should change to the correct answer.
  objective <- Maximize(c[1]*x[1] + c[2]*x[2])
  constraints <- list(x[1] <= h[1],
                      x[2] <= h[2],
                      A %*% x == b)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(result$value, 4)
  expect_equal(result$getValue(x), matrix(c(1, 2)), tolerance = TOL)
})

test_that("Test a basic LP with XPRESS", {
  skip_on_cran()
  skip_if_not(XPRESS_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
  result <- solve(prob, solver = "XPRESS")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "XPRESS")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  # XPRESS's default lower bound for a decision variable is zero
  # This quick test ensures that the cvxpy interface for XPRESS does *not* have that bound
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
})

test_that("Test a basic SOCP with XPRESS", {
  skip_if_not(XPRESS_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,2) + 1.0), list(x == 0))
  result <- solve(prob, solver = "XPRESS")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "XPRESS")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "XPRESS")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
})

test_that("Make sure XPRESS's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(XPRESS_INSTALLED)
  constraints <- list(x == 0)
  prob <- Problem(Minimize(p_norm(x,1)))
  result <- solve(prob, solver = "XPRESS")
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
  result <- solve(prob, solver = "XPRESS")
  duals_mosek <- lapply(constraints, function(c) { result$getDualValue(c) })
  result <- solve(prob, solver = "ECOS")
  duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
  for(i in seq_along(constraints))
    expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = TOL)
})

test_that("Test a basic LP with NAG", {
  skip_on_cran()
  skip_if_not(NAG_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,1) + 1.0), list(x == 0))
  result <- solve(prob, solver = "NAG")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "NAG")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  # NAG's default lower bound for a decision variable is zero
  # This quick test ensures that the cvxpy interface for NAG does *not* have that bound
  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "NAG")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
})

test_that("Test a basic SOCP with NAG", {
  skip_on_cran()
  skip_if_not(NAG_INSTALLED)
  prob <- Problem(Minimize(p_norm(x,2) + 1.0), list(x == 0))
  result <- solve(prob, solver = "NAG")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(0, 0)), tolerance = TOL)

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "NAG")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)

  objective <- Minimize(x[1])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "NAG")
  expect_equal(result$getValue(x), as.matrix(c(-100, 1)), tolerance = TOL)
})

test_that("Make sure NAG's dual result matches other solvers", {
  skip_on_cran()
  skip_if_not(NAG_INSTALLED)
  constraints <- list(x == 0)
  prob <- Problem(Minimize(p_norm(x,1)))
  result <- solve(prob, solver = "NAG")
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
  result <- solve(prob, solver = "NAG")
  duals_mosek <- lapply(constraints, function(c) { result$getDualValue(c) })
  result <- solve(prob, solver = "ECOS")
  duals_ecos <- lapply(constraints, function(c) { result$getDualValue(c) })
  for(i in seq_along(constraints))
    expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = TOL)
})

test_that("Test the list of installed solvers", {
  skip_on_cran()
  prob <- Problem(Minimize(p_norm(x, 1) + 1.0), list(x == 0))
  for(solver in names(SOLVER_MAP_CONIC)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$value, 1.0, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(0, 0)), tolerance = TOL)
    } else
      expect_error(result <- solve(prob, solver = solver), paste("The solver", solver, "is not installed"))
  }

  for(solver in names(SOLVER_MAP_QP)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$getValue(x), matrix(c(0, 0)), tolerance = TOL)
    } else
      expect_error(result <- solve(prob, solver = solver), paste("The solver", solver, "is not installed"))
  }
})
