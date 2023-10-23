context("test-g01-qp_solvers")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")
w <- Variable(5, name = "w")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

slope <- Variable(1, name = "slope")
offset <- Variable(1, name = "offset")
quadratic_coeff <- Variable(1, name = "quadratic_coeff")

Tlen <- 30
position <- Variable(2, Tlen, name = "position")
velocity <- Variable(2, Tlen, name = "velocity")
force <- Variable(2, Tlen - 1, name = "force")

xs <- Variable(80, name = "xs")
xsr <- Variable(50, name = "xsr")
xef <- Variable(80, name = "xef")

# Check for all installed QP solvers
solvers <- installed_solvers()
solvers <- solvers[solvers %in% CVXR:::QP_SOLVERS]
# solvers  <-  setdiff(solvers, "CPLEX")   # Skip CPLEX on CRAN, since I believe there are some false positive failures.

solve_QP <- function(problem, solver_name) {
  solve(problem, solver = solver_name, verbose = TRUE)
}

test_quad_over_lin <- function(solver) {
  skip_on_cran()
  p <- Problem(Minimize(0.5*quad_over_lin(abs(x-1), 1)), list(x <= -1))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), matrix(c(-1, -1)), tolerance = 1e-4)

  for(con in constraints(p))
    expect_equal(result$getDualValue(con), matrix(c(2, 2)), tolerance = 1e-4)
}

test_abs <- function(solver) {
  skip_on_cran()
  u <- Variable(2)
  constr <- list()
  constr <- c(constr, abs(u[2] - u[1]) <= 100)
  prob <- Problem(Minimize(sum_squares(u)), constr)
  print(paste("The problem is QP:", is_qp(prob)))
  expect_true(is_qp(prob))
  result <- solve(prob, solver = solver)
  expect_equal(result$value, 0, tolerance = TOL)
}

test_power <- function(solver) {
  skip_on_cran()
  p <- Problem(Minimize(sum(power(x, 2))), list())
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), matrix(c(0, 0)), tolerance = 1e-4)
}

test_power_matrix <- function(solver) {
  skip_on_cran()
  p <- Problem(Minimize(sum(power(A - 3, 2))), list())
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), matrix(3, nrow = 2, ncol = 2), tolerance = 1e-4)
}

test_square_affine <- function(solver) {
  skip_on_cran()
  A <- matrix(rnorm(10*2), nrow = 10, ncol = 2)
  b <- matrix(rnorm(10), nrow = 10, ncol = 1)
  p <- Problem(Minimize(sum_squares(A %*% x - b)))
  result <- solve_QP(p, solver)

  x_star <- base::solve(qr(A),b)
  for(var in variables(p))
    expect_equal(result$getValue(var), x_star, tolerance = 0.1)
}

test_quad_form <- function(solver) {
  skip_on_cran()
  set.seed(0)
  A <- matrix(rnorm(5*5), nrow = 5, ncol = 5)
  z <- matrix(rnorm(5), nrow = 5, ncol = 1)
  P <- t(A) %*% A
  q <- -2*P %*% z
  p <- Problem(Minimize(quad_form(w, P) + t(q) %*% w))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), z, tolerance = 1e-4)
}

test_rep_quad_form <- function(solver) {
  # A problem where the quad_form term is used multiple times.
  set.seed(0)
  A <- matrix(rnorm(5*5), nrow = 5, ncol = 5)
  z <- matrix(rnorm(5))
  P <- t(A) %*% A
  q <- -2*P %*% z
  qf <- QuadForm(w, P)
  p <- Problem(Minimize(0.5*qf + 0.5*qf + t(q) %*% w))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(z, result$getValue(var), tolerance = 1e-4)
}

test_affine_problem <- function(solver) {
  skip_on_cran()
  A <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  A <- pmax(A, 0)
  b <- matrix(rnorm(5), nrow = 5, ncol = 1)
  b <- pmax(b, 0)
  p <- Problem(Minimize(sum_entries(x)), list(x >= 0, A %*% x <= b))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), matrix(c(0, 0)), tolerance = 1e-3)
}

test_maximize_problem <- function(solver) {
  skip_on_cran()
  A <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  A <- pmax(A, 0)
  b <- matrix(rnorm(5), nrow = 5, ncol = 1)
  b <- pmax(b, 0)
  p <- Problem(Maximize(-sum_entries(x)), list(x >= 0, A %*% x <= b))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), matrix(c(0, 0)), tolerance = 1e-3)
}

test_norm_2 <- function(solver) {
  skip_on_cran()
  A <- matrix(rnorm(10*5), nrow = 10, ncol = 5)
  b <- matrix(rnorm(10), nrow = 10, ncol = 1)
  p <- Problem(Minimize(p_norm(A %*% w - b, 2)))
  result <- solve_QP(p, solver)

  x_star <- base::solve(qr(A),b)
  for(var in variables(p))
    expect_equal(result$getValue(var), x_star, tolerance = 0.1)
}

test_mat_norm_2 <- function(solver) {
  skip_on_cran()
  A <- matrix(rnorm(5*3), nrow = 5, ncol = 3)
  B <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  p <- Problem(Minimize(p_norm(A %*% C - B, 2)))
  result <- solve_QP(p, solver)

  C_star <- base::solve(qr(A), B)
  for(var in variables(p))
    expect_equal(result$getValue(var), C_star, tolerance = 0.1)
}

test_quad_form_coeff <- function(solver) {
  skip_on_cran()
  set.seed(0)
  A <- matrix(rnorm(5*5), nrow = 5, ncol = 5)
  z <- matrix(rnorm(5), nrow = 5, ncol = 1)
  P <- t(A) %*% A
  q <- -2*P %*% z
  p <- Problem(Minimize(quad_form(w, P) + t(q) %*% w))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), z, tolerance = 1e-4)
}

test_quad_form_bound <- function(solver) {
  skip_on_cran()
  P <- rbind(c(13, 12, -2), c(12, 17, 6), c(-2, 6, 12))
  q <- matrix(c(-22, -14.5, 13), nrow = 3, ncol = 1)
  r <- 1
  y_star <- matrix(c(1, 0.5, -1), nrow = 3, ncol = 1)
  p <- Problem(Minimize(0.5*quad_form(y, P) + t(q) %*% y + r),
               list(y >= -1, y <= 1))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), y_star, tolerance = 1e-4)
}

test_regression_1 <- function(solver) {
  skip_on_cran()
  set.seed(1)
  
  # Number of examples to use.
  n <- 100

  # Specify the true value of the variable.
  true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)

  # Generate data.
  x_data <- 5*rnorm(n)
  x_data_expanded <- sapply(1:3, function(i) { x_data^i })
  x_data_expanded <- matrix(x_data_expanded)
  y_data <- t(x_data_expanded) %*% true_coeffs + 0.5*runif(n)

  line <- offset + slope*x_data
  residuals <- line - y_data
  fit_error <- sum_squares(residuals)
  p <- Problem(Minimize(fit_error), list())
  result <- solve_QP(p, solver)

  model <- lm(y_data ~ x_data)
  expect_equal(sum((model$residuals)^2), result$value, tolerance = 1e-4)
  expect_equal(as.numeric(model$coefficients[1]), result$getValue(offset), tolerance = 1e-4)
  expect_equal(as.numeric(model$coefficients[2]), result$getValue(slope), tolerance = 1e-4)
}

test_regression_2 <- function(solver) {
  skip_on_cran()
  set.seed(1)
  
  # Number of examples to use.
  n <- 100

  # Specify the true value of the variable.
  true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)

  # Generate data.
  x_data <- 5*rnorm(n)
  x_data_expanded <- sapply(1:3, function(i) { x_data^i })
  x_data_expanded <- matrix(x_data_expanded)
  print(dim(x_data_expanded))
  print(dim(true_coeffs))
  y_data <- t(x_data_expanded) %*% true_coeffs + 0.5*runif(n)

  quadratic <- offset + slope*x_data + quadratic_coeff*x_data^2
  residuals <- quadratic - y_data
  fit_error <- sum_squares(residuals)
  p <- Problem(Minimize(fit_error), list())
  result <- solve_QP(p, solver)

  x_data_sq <- x_data^2
  model <- lm(y_data ~ x_data + x_data_sq)
  expect_equal(sum((model$residuals)^2), result$value, tolerance = 1e-4)
  expect_equal(as.numeric(model$coefficients[1]), result$getValue(offset), tolerance = 1e-4)
  expect_equal(as.numeric(model$coefficients[2]), result$getValue(slope), tolerance = 1e-4)
  expect_equal(as.numeric(model$coefficients[3]), result$getValue(quadratic_coeff), tolerance = 1e-4)
}

test_control <- function(solver) {
  skip_on_cran()
  
  # Some constraints on our motion.
  # The object should start from the origin, and end at rest.
  initial_velocity <- c(-20, 100)
  final_position <- c(100, 100)

  Tlen <- 30       # The number of timesteps.
  h <- 0.1         # The time between time intervals.
  mass <- 1        # Mass of object.
  drag <- 0.1      # Drag on object.
  g <- c(0, -9.8)  # Gravity on object.

  # Create a problem instance.
  constraints = list()

  # Add constraints on our variables.
  for(i in 1:(Tlen - 1)) {
    constraints <- c(constraints, position[,i+1] == position[,i] + h*velocity[,i])
    acceleration <- force[,i]/mass + g - drag*velocity[,i]
    constraints <- c(constraints, velocity[,i+1] == velocity[,i] + h*acceleration)
  }

  # Add position constraints.
  constraints <- c(constraints, position[,1] == 0)
  constraints <- c(constraints, position[,ncol(position)] == final_position)

  # Add velocity constraints.
  constraints <- c(constraints, velocity[,1] == initial_velocity)
  constraints <- c(constraints, velocity[,ncol(velocity)] == 0)

  # Solve the problem.
  p <- Problem(Minimize(0.01*sum_squares(force)), constraints)
  result <- solve_QP(p, solver)
  expect_equal(1059.616, result$value, tolerance = 0.1)
}

test_sparse_system <- function(solver) {
  skip_on_cran()
  if(!require(Matrix)) {
    print("Matrix library not installed. Skipping test.")
    return()
  }
  
  m <- 100
  n <- 80
  
  set.seed(1)
  density <- 0.4
  A <- rsparsematrix(m, n, density)
  b <- rnorm(m)

  p <- Problem(Minimize(sum_squares(A %*% xs - b)), list(xs == 0))
  result <- solve_QP(p, solver)
  expect_equal(sum(b^2), result$value, tolerance = 1e-4)
}

test_smooth_ridge <- function(solver) {
  skip_on_cran()
  set.seed(1)
  n <- 50
  k <- 20
  eta <- 1

  A <- matrix(1, nrow = k, ncol = n)
  b <- rep(1, k)
  obj <- sum_squares(A %*% xsr - b) + eta*sum_squares(diff(xsr))
  p <- Problem(Minimize(obj), list())
  result <- solve_QP(p, solver)
  expect_equal(result$value, 0, tolerance = 1e-4)
}

test_huber_small <- function(solver) {
  skip_on_cran()
  
  # Solve the Huber regression problem.
  x <- Variable(3)
  objective <- sum(huber(x))

  # Solve problem with QP.
  p <- Problem(Minimize(objective), list(x[3] >= 3))
  result <- solve_QP(p, solver)
  expect_equal(result$getValue(x[3]), 3, tolerance = 1e-4)
  expect_equal(result$getValue(objective), 5, tolerance = 1e-4)
}

test_huber <- function(solver) {
  skip_on_cran()
  if(!require(Matrix)) {
    print("Matrix library not installed. Skipping test.")
    return()
  }

  # Generate problem data.
  n <- 3
  m <- 5

  set.seed(2)
  A <- rsparsematrix(m, n, density = 0.8)
  x_true <- matrix(rnorm(n), ncol = 1)/sqrt(n)
  ind95 <- as.numeric(rnorm(m) < 0.95)
  b <- A %*% x_true + 0.5*rnorm(m)*ind95 + 10*runif(m)*(1 - ind95)

  # Solve the Huber regression problem.
  x <- Variable(n)
  objective <- sum(CVXR::huber(A %*% x - b))

  # Solve problem with QP.
  p <- Problem(Minimize(objective))
  result <- solve_QP(p, solver)
  expect_equal(result$getValue(objective), 1.327429461061672, tolerance = 1e-3)
  expect_equal(result$getValue(x), matrix(c(-1.03751745, 0.86657204, -0.9649172)), tolerance = 1e-3)
}

test_equivalent_forms_1 <- function(solver) {
  skip_on_cran()
  m <- 100
  n <- 80
  r <- 70

  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  G <- matrix(rnorm(r*n), nrow = r, ncol = n)
  h <- matrix(rnorm(r), nrow = r, ncol = 1)

  obj1 <- 0.1*sum((A %*% xef - b)^2)
  cons <- list(G %*% xef == h)

  p1 <- Problem(Minimize(obj1), cons)
  result <- solve_QP(p1, solver)
  expect_equal(result$value, 68.1119420108, tolerance = 1e-4)
}

test_equivalent_forms_2 <- function(solver) {
  skip_on_cran()
  m <- 100
  n <- 80
  r <- 70

  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  G <- matrix(rnorm(r*n), nrow = r, ncol = n)
  h <- matrix(rnorm(r), nrow = r, ncol = 1)

  # ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2.
  P <- t(A) %*% A
  q <- -2*t(A) %*% b
  r <- t(b) %*% b

  obj2 <- 0.1*(quad_form(xef, P) + t(q) %*% xef + r)
  cons <- list(G %*% xef == h)

  p2 <- Problem(Minimize(obj2), cons)
  result <- solve_QP(p2, solver)
  expect_equal(result$value, 68.1119420108, tolerance = 1e-4)
}

test_equivalent_forms_3 <- function(solver) {
  skip_on_cran()
  m <- 100
  n <- 80
  r <- 70

  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  G <- matrix(rnorm(r*n), nrow = r, ncol = n)
  h <- matrix(rnorm(r), nrow = r, ncol = 1)

  # ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2.
  P <- t(A) %*% A
  q <- -2*t(A) %*% b
  r <- t(b) %*% b
  Pinv <- base::solve(P)

  obj3 <- 0.1 * (matrix_frac(xef, Pinv) + t(q) %*% xef + r)
  cons <- list(G %*% xef == h)

  p3 <- Problem(Minimize(obj3), cons)
  result <- solve_QP(p3, solver)
  expect_equal(result$value, 68.1119420108, tolerance = 1e-4)
}

test_that("test all solvers", {
  skip_on_cran()
  for(solver in solvers) {
    test_quad_over_lin(solver)
    test_power(solver)
    test_power_matrix(solver)
    test_square_affine(solver)
    test_quad_form(solver)
    test_affine_problem(solver)
    test_maximize_problem(solver)
    test_abs(solver)

    # Do we need the following functionality?
    # test_norm_2(solver)
    # test_mat_norm_2(solver)

    test_quad_form_coeff(solver)
    test_quad_form_bound(solver)
    test_regression_1(solver)
    test_regression_2(solver)
    test_rep_quad_form(solver)

    # Slow tests:
    test_control(solver)
    test_sparse_system(solver)
    test_smooth_ridge(solver)
    test_huber_small(solver)
    test_huber(solver)
    test_equivalent_forms_1(solver)
    test_equivalent_forms_2(solver)
    test_equivalent_forms_3(solver)
  }
})

test_that("Test warm start", {
  skip_on_cran()
  m <- 200
  n <- 100

  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- Parameter(m)

  # Construct the problem.
  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A %*% x - b)))

  value(b) <- rnorm(m)
  result <- solve(prob, solver = "OSQP", warm_start=FALSE)
  result2 <- solve(prob, solver = "OSQP", warm_start=TRUE)
  expect_equal(result$value, result2$value)

  value(b) <- rnorm(m)
  result <- solve(prob, solver = "OSQP", warm_start=TRUE)
  result2 <- solve(prob, solver = "OSQP", warm_start=FALSE)
  expect_equal(result$value, result2$value)
})

test_that("Test GUROBI warm start with a user provided point", {
  skip_on_cran()
  if("GUROBI" %in% installed_solvers()) {
    library(gurobi)
    m <- 4
    n <- 3
    
    y <- Variable(nonneg = TRUE)
    X <- Variable(m, n)
    X_vals <- matrix(0:(m*n - 1), nrow = m, ncol = n, byrow = TRUE)
    value(X) <- X_vals + 1
    prob <- Problem(Minimize(y^2 + sum(x)), list(X == X_vals))
    result <- solve(prob, solver = "GUROBI", warm_start = TRUE)
    
    # Check that "start" value was set appropriately.
    model <- result$solver_stats@extra_stats
    model_x <- getVars(model)
    expect_equal(model_x[1]$start, gurobi:::UNDEFINED)
    expect_true(is.allclose(model_x[1]$x, 0))
    for(i in 1:size(X)) {
      row <- (i - 1) %% nrow(X) + 1
      col <- floor((i - 1) / nrow(X)) + 1
      expect_equal(model_x[i + 1]$start, X_vals[row, col] + 1)
      expect_true(is.allclose(model_x[i + 1]$x, result$getValue(X)[row, col]))
    }
  }
})

test_that("Test solve parametric vs. full problem", {
  skip_on_cran()
  x <- Variable()
  a <- 10
  # b_vec <- c(-10, -2., 2., 3., 10.)
  b_vec <- c(-10, -2.)

  for(solver in solvers) {
    # Solve from scratch with no parameters.
    x_full <- list()
    obj_full <- c()
    for(b in b_vec) {
      obj <- Minimize(a*x^2 + b*x)
      constraints <- list(0 <= x, x <= 1)
      prob <- Problem(obj, constraints)
      result <- solve(prob, solver=solver)
      x_full <- c(x_full, list(result$getValue(x)))
      obj_full <- c(obj_full, result$value)
    }

    # Solve parametric.
    x_param <- list()
    obj_param <- c()
    b <- Parameter()
    obj <- Minimize(a*x^2 + b*x)
    constraints <- list(0 <= x, x <= 1)
    prob <- Problem(obj, constraints)
    for(b_value in b_vec) {
      value(b) <- b_value
      result <- solve(prob, solver=solver)
      x_param <- c(x_param, list(result$getValue(x)))
      obj_param <- c(obj_param, result$value)
    }

    print(x_full)
    print(x_param)
    for(i in seq_along(b_vec)) {
      expect_equal(x_full[[i]], x_param[[i]], tolerance = 1e-3)
      expect_equal(obj_full[i], obj_param[i], tolerance = TOL)
    }
  }
})

test_that("Test issue arising with square plus parameter", {
  skip_on_cran()
  a <- Parameter(value=1)
  b <- Variable()

  obj <- Minimize(b^2 + abs(a))
  prob <- Problem(obj)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 1.0, tolerance = TOL)
})

# TODO: Update tests for the R gurobi library.

test_that("Test GUROBI time limit no solution", {
  # Make sure that if Gurobi terminates due to a time limit before finding a solution:
  #   1) no error is raised,
  #   2) solver stats are returned.
  # The test is skipped if something changes on Gurobi's side so that:
  #   - a solution is found despite a time limit of zero,
  #   - a different termination criteria is hit first.
  skip_on_cran()
  
  if("GUROBI" %in% installed_solvers()) {
    library(gurobi)
    objective <- Minimize(x[1])
    constraints <- list(x[1] >= 1)
    prob <- Problem(objective, constraints)
    tryCatch({
      result <- solve(prob, solver = "GUROBI", TimeLimit = 0.0)
    }, error = function(e) {
      stop("An exception ", e, " is raised instead of returning a result.")
    })
    
    extra_stats <- NULL
    solver_stats <- ifelse("solver_stats" %in% names(result), result$solver_stats, NULL)
    if(!is.null(solver_stats))
      extra_stats <- ifelse("extra_stats" %in% slotNames(solver_stats), solver_stats@extra_stats, NULL)
    expect_true(!is.null(extra_stats), "Solver stats have not been returned")
    
    nb_solutions <- ifelse("SolCount" %in% slotNames(extra_stats), extra_stats@SolCount, NULL)
    if(!is.null(nb_solutions)) {
      print("GUROBI has found a solution, the test is not relevant anymore")
      return()
    }
    
    solver_status <- ifelse("Status" %in% slotNames(extra_stats), extra_stats@Status, NA_character_)
    if(solver_status != gurobi:::TIME_LIMIT) {
      print("GUROBI terminated for a different reason than reaching time limit, the test is not relevant anymore")
      return()
    }
  } else {
    prob <- Problem(Minimize(norm1(x)), list(x == 0))
    expect_error(solve(prob, solver = "GUROBI", TimeLimit = 0), "The solver GUROBI is not installed.", fixed = TRUE)
  }
})

test_that("test GUROBI environment", {
  # Tests that Gurobi environments can be passed to Model.
  # Gurobi environments can include licensing and model parameter data.
  skip_on_cran()
  
  if("GUROBI" %in% installed_solvers()) {
    library(gurobi)
    
    # Set a few parameters to random values close to their defaults.
    params <- list(MIPGap = runif(),                  # range {0, INFINITY}
                   AggFill = sample(0:9, size = 1),   # range {-1, MAXINT}
                   PerturbValue = runif())            # range: {0, INFINITY}
    
    # Create a custom environment and set some parameters.
    custon_env <- gurobi:::Env()
    for(k in names(params)) {
      v <- params[[k]]
      custom_env$setParam(k, v)
    }
    
    # Testing QP Solver interface.
    sth <- StandardTestLPs.test_lp_0(solver = "GUROBI", env = custom_env)
    model <- sth@prob@solver_stats@extra_stats
    for(k in names(params)) {
      # https://www.gurobi.com/documentation/current/refman/cpp_env_getparaminfo.html
      v <- params[[k]]
      info <- model$getParamInfo(k)
      p_val <- info[[3]]
      expect_equal(v, p_val)
    }
  } else {
    prob <- Problem(Minimize(norm1(x)), list(x == 0))
    expect_error(solve(prob, solver = "GUROBI", TimeLimit = 0), "The solver GUROBI is not installed.", fixed = TRUE)
  }
})
