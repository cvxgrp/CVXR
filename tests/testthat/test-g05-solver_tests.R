INSTALLED_SOLVERS <- installed_solvers()
TOL <- 1e-4
x <- Variable(2, name = "x")
slope <- Variable(1, name = "slope")
offset <- Variable(1, name = "offset")
quadratic_coeff <- Variable(1, name = "quadratic_coeff")

#LP, QP, SOCP, SDP, MIP, GP
test_that("Test a basic LP with all solvers", {
  skip_on_cran()

  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, x[1] + 2*x[2] <= 3, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)

  for(solver in INSTALLED_SOLVERS){
    result <- solve(prob, solver = solver)
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
  }

  #Example from
  # https://www.cvxpy.org/examples/basic/linear_program.html
  m <- 15
  n <- 10
  set.seed(1)
  s0 <- rnorm(m)
  lamb0 <- pmax(-s0, 0)
  s0 <- pmax(s0, 0)
  x0 <- rnorm(n)
  A <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b <- A %*% x0 + s0
  c <- -t(A) %*% lamb0
  x_var <- Variable(n)
  prob <- Problem(Minimize(t(c) %*% x_var), list(A %*% x_var <= b))
  for(solver in INSTALLED_SOLVERS){
    suppressWarnings(
    result <- solve(prob, solver = solver, abstol = 1e-6)
    )
    expect_equal(result$value, -7.83446, tolerance = TOL)
  }

})

#QP
test_that("Test a basic QP with all solvers", {
  skip_on_cran()

  set.seed(0)
  #Regression 1 test-g01-qp.R
  n <- 100

  # Specify the true value of the variable
  true_coeffs <- matrix(c(2, -2, 0.5), ncol = 1)

  # Generate data
  x_data <- 5*rnorm(n)
  x_data_expanded <- sapply(1:3, function(i) { x_data^i })
  y_data <- x_data_expanded %*% true_coeffs + 0.5*runif(n)

  line <- offset + slope*x_data
  residuals <- line - y_data
  fit_error <- sum_squares(residuals)
  p <- Problem(Minimize(fit_error), list())
  model <- lm(y_data ~ x_data)

  ##CBC, ECOS_BB, GLPK_MI, GLPK don't support QPs
  applicable_solvers  <- setdiff(INSTALLED_SOLVERS, list("CBC", "ECOS_BB", "GLPK_MI", "GLPK"))
  ## SCS is known to cause problems
  ## Works sometimes if we set rho_x = 1e-2, but too problematic
  applicable_solvers  <- setdiff(applicable_solvers, c("SCS", "CVXOPT"))

  for(solver in applicable_solvers) {
    result <- solve(p, solver)
    print(sprintf("Solver: %s status: %s\n", solver, result$status))
    expect_equal(sum((model$residuals)^2), result$value, tolerance = TOL)
    expect_equal(as.numeric(model$coefficients[1]), result$getValue(offset), tolerance = TOL)
    expect_equal(as.numeric(model$coefficients[2]), result$getValue(slope), tolerance = TOL)
  }

  ##Regression 2 test-g01-qp.R
  quadratic <- offset + slope*x_data + quadratic_coeff*x_data^2
  residuals <- quadratic - y_data
  fit_error <- sum_squares(residuals)
  p <- Problem(Minimize(fit_error), list())
  x_data_sq <- x_data^2
  model <- lm(y_data ~ x_data + x_data_sq)
  for(solver in applicable_solvers) {
    ##cat(sprintf("Doing %s\n", solver))
    result <- solve(p, solver)
    expect_equal(sum((model$residuals)^2), result$value, tolerance = TOL)
    expect_equal(as.numeric(model$coefficients[1]), result$getValue(offset), tolerance = TOL)
    expect_equal(as.numeric(model$coefficients[2]), result$getValue(slope), tolerance = TOL)
    expect_equal(as.numeric(model$coefficients[3]), result$getValue(quadratic_coeff), tolerance = TOL)
  }

})

test_that("Test a basic SOCP with solvers", {
  skip_on_cran()
  applicable_solvers <- setdiff(INSTALLED_SOLVERS, list("CBC", "ECOS_BB", "GLPK_MI", "GLPK", "OSQP"))
  # Example from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3, (x[1] + 2*x[2])^2 <= 9, x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)
  for(solver in applicable_solvers) {
    result <- solve(prob, solver = solver)
    print(sprintf("Solver: %s status: %s\n", solver, result$status))
    expect_equal(result$value, -9, tolerance = TOL)
    expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
  }

  # Example from
  # https://www.cvxpy.org/examples/basic/socp.html
  set.seed(1)
  m <- 3
  n <- 10
  p <- 5
  n_i <- 5
  f = rnorm(n)
  A <- vector(mode="list", m)
  b <- vector(mode="list", m)
  c <- vector(mode="list", m)
  d <- vector(mode="list", m)
  x0 <- rnorm(n)
  soc_constraints <- vector(mode="list", m+1)
  x_var <- Variable(n)
  for(i in 1:m){
    A[[i]] <- matrix(rnorm(n_i * n), nrow = n_i, ncol = n)
    b[[i]] <- rnorm(n_i)
    c[[i]] <- rnorm(n)
    d[[i]] <- norm(A[[i]] %*% x0 + b[[i]], "f")
    soc_constraints[[i]] <- CVXR:::SOC(t(c[[i]]) %*% x_var + d[[i]], A[[i]] %*% x_var + b[[i]])
  }
  Fmat <- matrix(rnorm(p * n), nrow = p, ncol = n)
  g <- Fmat %*% x0
  soc_constraints[[m+1]] <- Fmat %*% x_var == g
  prob <- Problem(Minimize(t(f)%*%x_var), soc_constraints)
  for(solver in applicable_solvers) {
    result <- solve(prob, solver=solver)
    print(sprintf("Solver: %s status: %s\n", solver, result$status))
    expect_equal(result$value, -10.13341, tolerance=TOL)
  }

})

test_that("Test a basic SDP with all solvers", {
  skip_on_cran()
  # Example from
  # https://www.cvxpy.org/examples/basic/sdp.html
  n <- 3
  p <- 3
  set.seed(4)
  C <- matrix(rnorm(n^2), nrow = n)
  A <- vector(mode="list", p)
  b <- numeric(p)
  X <- Variable(n, n, PSD=T)
  constraints <- list()
  for(i in 1:p){
    A[[i]] <- matrix(rnorm(n^2), nrow = n)
    b[i] <- rnorm(1)
    constraints <- c(constraints, matrix_trace(A[[i]] %*% X) == b[i])
  }
  prob <- Problem(Minimize(matrix_trace(C %*% X)), constraints)

  check_matrix <- matrix(c(.3771707, -.5692205, .1166577, -.5692205,
                           .8590441, -.1760618, .11665767,
                           -.17606183, .03608155), nrow = n)
  # SCS and MOSEK and CVXOPT are the only three solvers that support SDPs
  for(solver in intersect(INSTALLED_SOLVERS, c("SCS", "MOSEK", "CVXOPT"))) {
    result <- solve(prob, solver = solver)
    print(sprintf("Solver: %s status: %s\n", solver, result$status))
    expect_equal(result$value, 1.395479, tolerance=TOL)
    expect_equal(result$getValue(X), check_matrix, tolerance=TOL)
  }
})

test_that("Test a mixed-integer quadratic program",{
  skip_on_cran()
  # Example from
  # https://www.cvxpy.org/examples/basic/mixed_integer_quadratic_program.html
  m <- 40
  n <- 25
  set.seed(1)
  A <- matrix(runif(m*n), nrow = m)
  b <- rnorm(m)
  x <- Variable(n, integer=T)
  obj <- Minimize(sum_squares(A %*% x - b))
  prob <- Problem(obj)
  applicable_solvers <- intersect(INSTALLED_SOLVERS, c("ECOS_BB", "CPLEX", "GUROBI", "MOSEK"))
  for(solver in applicable_solvers) {
    result <- solve(prob, solver=solver)
    expect_equal(result$value, 13.34086, tolerance=TOL)
  }
})

test_that("Test a simple geometric program", {
  skip_on_cran()
  # Example from
  # https://www.cvxpy.org/examples/dgp/dgp_fundamentals.html
  x <- Variable(pos=TRUE)
  y <- Variable(pos=TRUE)
  z <- Variable(pos=TRUE)

  obj <- x * y * z
  constraints <- list(4*x*y*z + 2*x*z <= 10,
                      x <= 2*y,
                      y <= 2*x,
                      z >=1)
  prob <- Problem(Maximize(obj), constraints)
  for(solver in c("ECOS", "SCS")){
    result <- solve(prob, solver=solver, gp=TRUE)
    expect_equal(result$value, 2, tolerance=TOL)
    expect_equal(result$getValue(x), 1, tolerance=TOL)
    expect_equal(result$getValue(y), 2, tolerance=TOL)
    expect_equal(result$getValue(z), 1, tolerance=TOL)
  }

})
