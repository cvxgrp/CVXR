TOL <- 1e-6

test_that("Find the largest Euclidean ball in the polyhedron", {
  # Create the input data
  a1 <- matrix(c(2,1), nrow = 2, ncol = 1)
  a2 <- matrix(c(2,-1), nrow = 2, ncol = 1)
  a3 <- matrix(c(-1,2), nrow = 2, ncol = 1)
  a4 <- matrix(c(-1,-2), nrow = 2, ncol = 1)
  b <- rep(1,4)
  
  # Create and solve the model
  r <- Variable(name = "r")
  x_c <- Variable(2, name = "x_c")
  obj <- Maximize(r)
  constraints <- list(
    t(a1) %*% x_c + norm(a1,"F")*r <= b[1],
    t(a2) %*% x_c + norm(a2,"F")*r <= b[2],
    t(a3) %*% x_c + norm(a3,"F")*r <= b[3],
    t(a4) %*% x_c + norm(a4,"F")*r <= b[4]
  )
  
  p <- Problem(obj, constraints)
  result <- solve(p)
  expect_equal(result$value, 0.447214, tolerance = TOL)
  expect_equal(result$getValue(r), result$value, tolerance = TOL)
  expect_equal(result$getValue(x_c), matrix(c(0,0)), tolerance = TOL)
})

test_that("Test examples from the README", {
  # Problem data
  m <- 30
  n <- 20
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  
  # Construct the problem
  x <- Variable(n)
  objective <- Minimize(sum((A %*% x - b)^2))
  constraints <- list(x >= 0, x <= 1)
  p <- Problem(objective, constraints)
  
  # The optimal objective is returned by solve(p)
  result <- solve(p)
  # The optimal value for x is stored in result$getValue(x)
  print(result$getValue(x))
  # The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraints[[1]])
  print(result$getDualValue(constraints[[1]]))

  ###########################################
  # Scalar variable
  a <- Variable()
  
  # Column vector variable of length 5
  x <- Variable(5)
  
  # Matrix variable with 4 rows and 7 columns
  A <- Variable(4, 7)
  
  ###########################################
  # Positive scalar parameter
  m <- Parameter(sign = "positive")
  
  # Column vector parameter with unknown sign (by default)
  c <- Parameter(5)
  
  # Matrix parameter with negative entries
  G <- Parameter(4, 7, sign = "negative")
  
  # Assigns a constant value to G
  value(G) <- -matrix(1, nrow = 4, ncol = 7)
  
  expect_error(value(G) <- matrix(1, nrow = 4, ncol = 7))
  
  ###########################################
  a <- Variable()
  x <- Variable(5)
  
  # expr is an Expression object after each assignment
  expr <- 2*x
  expr <- expr - a
  expr <- sum(expr) + norm2(x)
  
  ###########################################
  # Problem data
  n <- 10
  m <- 5
  A <- matrix(rnorm(n*m), nrow = n, ncol = m)
  b <- matrix(rnorm(n), nrow = n, ncol = 1)
  gamma <- Parameter(sign = "positive")
  
  # Construct the problem
  x <- Variable(m)
  loss <- sum((A %*% x - b)^2)
  
  get_x <- function(gamma) {
    objective <- Minimize(loss + gamma*norm1(x))
    p <- Problem(objective)
    
    result <- solve(p)
    result$getValue(x)
  }
  
  gammas <- 10^seq(-1, 2, length.out = 2)
  
  # Serial computation
  x_values <- sapply(gammas, get_x)
  
  ###########################################
  n <- 10
  
  mu <- matrix(rnorm(n), nrow = 1, ncol = n)
  sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
  sigma <- t(sigma) %*% sigma
  gamma <- Parameter(sign = "positive")
  value(gamma) <- 1
  x <- Variable(n)
  
  # Constants:
  # mu is the vector of expected returns.
  # sigma is the covariance matrix.
  # gamma is a Parameter that trades off risk and return.
  
  # Variables:
  # x i s a vector of stock holdings as fractions of total assets.
  
  expected_return <- mu %*% x
  risk <- quad_form(x, sigma)
  
  gamma <- 1
  objective <- Maximize(expected_return - gamma*risk)
  p <- Problem(objective, list(sum(x) == 1))
  result <- solve(p)
  
  # The optimal expected return
  print(result$getValue(expected_return))
  
  # The optimal risk
  print(result$getValue(risk))
  
  ###########################################
  N <- 50
  M <- 40
  n <- 10
  data1 <- lapply(1:N, function(i) { list(1, matrix(rnorm(n, mean = 1.0, sd = 2.0))) })
  data2 <- lapply(1:M, function(i) { list(-1, matrix(rnorm(n, mean = -1.0, sd = 2.0))) })
  data <- c(data1, data2)
  
  # Construct problem
  gamma <- Parameter(sign = "positive")
  value(gamma) <- 0.1
  a <- Variable(n)
  b <- Variable()
  
  slack <- lapply(data, function(x) { 
    label <- x[[1]]
    sample <- x[[2]]
    Pos(1 - label*(t(sample) %*% a - b))
  })
  gamma <- 0.1
  objective <- Minimize(norm2(a) + gamma*Reduce("+", slack))
  p <- Problem(objective)
  result <- solve(p)
  
  # Count misclassifications
  errors <- 0
  for(v in data) {
    label <- v[[1]]
    sample <- v[[2]]
    if(label * result$getValue(t(sample) %*% a - b) < 0)
      errors <- errors + 1
  }
  
  print(errors)
  print(result$getValue(a))
  print(result$getValue(b))
})

test_that("Test advanced tutorial", {
  # Solving a problem with different solvers
  x <- Variable(2)
  obj <- Minimize(x[1] + norm1(x))
  constraints <- list(x >= 2)
  prob <- Problem(obj, constraints)
  
  # Solve with ECOS
  result <- solve(prob, solver = "ECOS")
  print(paste("optimal value with ECOS:", result$value))
  expect_equal(result$value, 6, tolerance = TOL)
  
  # Solve with SCS
  result <- solve(prob, solver = "SCS")
  print(paste("optimal value with SCS:", result$value))
  expect_equal(result$value, 6, tolerance = 1e-3)
  
  x <- Variable()
  prob <- Problem(Minimize(x^2), list(x == 2))
  
  # Get ECOS arguments
  data <- get_problem_data(prob, "ECOS")
  
  # Get SCS arguments
  data <- get_problem_data(prob, "SCS")
})

test_that("Test log determinant", {
  x <- t(data.frame(c(0.55, 0.25, -0.2, -0.25, -0.0, 0.4),
                  c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)))
  n <- nrow(x)
  m <- ncol(x)
  
  # Create and solve the model
  A <- Variable(n, n)
  b <- Variable(n)
  obj <- Maximize(log_det(A))
  constraints <- lapply(1:m, function(i) { p_norm(A %*% as.matrix(x[,i]) + b) <= 1 })
  p <- Problem(obj, constraints)
  result <- solve(p)
  expect_equal(result$value, 1.9746, tolerance = 1e-2)
})

test_that("Test portfolio problem", {
  require(Matrix)
  set.seed(5)
  n <- 100
  m <- 10
  pbar <- matrix(1, nrow = n, ncol = 1) * 0.03 +
          c(matrix(rnorm(n-1), nrow = n-1, ncol = 1), 0) * 0.12
  
  Fmat <- rsparsematrix(m, n, density = 0.01)
  Fmat@x <- rep(1, length(Fmat@x))
  D <- sparseMatrix(1:n, 1:n, x = rep(1,n))
  D@x <- rnorm(length(D@x))^2
  Z <- matrix(rnorm(m), nrow = m, ncol = 1)
  Z <- Z %*% t(Z)
  
  x <- Variable(n)
  y <- Fmat %*% x
  mu <- 1
  ret <- t(pbar) %*% x
  # DCP attr causes error because not all the curvature
  # matrices are reduced to constants when an atom is scalar.
  risk <- p_norm(D %*% x)^2 + (Z %*% y)^2
})

test_that("Test examples from CVXR introduction", {
  m <- 30
  n <- 20
  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  
  # Construct the problem
  x <- Variable(n)
  objective <- Minimize(sum_squares(A %*% x - b))
  constraints <- list(0 <= x, x <= 1)
  prob <- Problem(objective, constraints)
  
  # The optimal objective is returned by solve(prob)
  result <- solve(prob)
  # The optimal value for x
  print(result$getValue(x))
  # The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraints[[1]])
  print(result$getDualValue(constraints[[1]]))
  
  ###########################################
  # Create two scalar variables
  x <- Variable()
  y <- Variable()
  
  # Create two constraints
  constraints <- list(x + y == 1, x - y >= 1)
  
  # Form objective
  obj <- Minimize((x - y)^2)
  
  # Form and solve problem
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  cat("status:", result$status)
  cat("\noptimal value:", result$value) 
  cat("\noptimal var:", result$getValue(x), ", ", result$getValue(y))
  
  expect_equal(tolower(result$status), "optimal")
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 0, tolerance = TOL)
  
  ###########################################
  # Replace the objective
  prob@objective <- Maximize(x + y)
  result <- solve(prob)
  cat("optimal value:", result$value)
  
  expect_equal(result$value, 1.0, tolerance = TOL)
  
  # Replace the constraints (x + y == 1)
  prob@constraints[[1]] <- (x + y <= 3)
  result <- solve(prob)
  cat("optimal value:", result$value)
  
  expect_equal(result$value, 3.0, tolerance = TOL)
  
  ###########################################
  x <- Variable()
  
  # An infeasible problem
  prob <- Problem(Minimize(x), list(x >= 1, x <= 0))
  result <- solve(prob)
  cat("status:", result$status)
  cat("optimal value:", result$value)
  
  expect_equal(tolower(result$status), "infeasible")
  expect_equal(result$value, Inf)
  
  # An unbounded problem
  prob <- Problem(Minimize(x))
  result <- solve(prob)
  cat("status:", result$status)
  cat("optimal value:", result$value)
  
  expect_equal(tolower(result$status), "unbounded")
  expect_equal(result$value, -Inf)
  
  ###########################################
  # A scalar variable
  a <- Variable()
  
  # Column vector variable of length 5
  x <- Variable(5)
  
  # Matrix variable with 4 rows and 7 columns
  A <- Variable(4, 7)
  
  ###########################################
  m <- 10
  n <- 5
  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  
  # Construct the problem
  x <- Variable(n)
  objective <- Minimize(sum((A %*% x - b)^2))
  constraints <- list(0 <= x, x <= 1)
  prob <- Problem(objective, constraints)
  
  result <- solve(prob)
  cat("Optimal value:", result$value)
  cat("Optimal var:", result$getValue(x))
  
  expect_equal(result$value, 7.244277, tolerance = TOL)
  
  ###########################################
  # Positive scalar parameter
  m <- Parameter(sign = "positive")
  
  # Column vector parameter with unknown sign (by default)
  c <- Parameter(5)
  
  # Matrix parameter with negative entries
  G <- Parameter(4, 7, sign = "negative")
  
  # Assigns a constant value to G
  value(G) <- -matrix(1, nrow = 4, ncol = 7)
  
  ###########################################
  # Create parameter, then assign value
  rho <- Parameter(sign = "positive")
  value(rho) <- 2
  
  # Initialize parameter with a value
  rho <- Parameter(sign = "positive", value = 2)
  
  ###########################################
  # Problem data
  n <- 15
  m <- 10
  set.seed(1)
  A <- matrix(rnorm(n*m), nrow = n, ncol = m)
  b <- matrix(rnorm(n), nrow = n, ncol = 1)
  # gamma must be positive due to DCP rules
  gamma <- Parameter(sign = "positive")
  
  # Construct the problem
  x <- Variable(m)
  error <- sum_squares(A %*% x - b)

  # Construct a trade-off curve of ||Ax-b||^2 vs. ||x||_1
  sq_penalty <- c()
  l1_penalty <- c()
  x_values <- list()
  gammas <- 10^seq(-4, 6, length.out = 50)
  for(gamma in gammas) {
    obj <- Minimize(error + gamma*p_norm(x,1))
    prob <- Problem(obj)

    result <- solve(prob)
    sq_penalty <- c(sq_penalty, result$getValue(error))
    l1_penalty <- c(l1_penalty, result$getValue(norm1(x)))
    x_values <- c(x_values, result$getValue(x))
  }
  
  ###########################################
  X <- Variable(5, 4)
  A <- matrix(1, nrow = 3, ncol = 5)
  
  # Use size(expr) to get the dimensions
  cat("dimensions of X:", size(X))
  cat("\ndimensions of sum_entries(X):", size(sum_entries(X)))
  cat("\ndimensions of A %*% X:", size(A %*% X))
  
  # ValueError raised for invalid dimensions
  expect_error(A + X)
})

test_that("Test image in-painting", {
  set.seed(1)
  rows <- 100
  cols <- 100
  
  # Load the images and convert to arrays
  Uorig <- matrix(sample(0:255, size = rows * cols, replace = TRUE), nrow = rows, ncol = cols)
  rows <- nrow(Uorig)
  cols <- ncol(Uorig)
  
  # Known is 1 if the pixel is known, 0 if the pixel was corrupted
  Known <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:rows) {
    for(j in 1:cols) {
      if(runif(1) > 0.7)
        Known[i,j] <- 1
    }
  }
  Ucorr <- Known %*% Uorig
  
  # Recover the original image using total variation in-painting
  U <- Variable(rows, cols)
  obj <- Minimize(tv(U))
  constraints <- list(Known * U == Known * Ucorr)
  prob <- Problem(obj, constraints)
  solve(prob, solver = "SCS")
})

test_that("Test the log_sum_exp function", {
  set.seed(1)
  m <- 5
  n <- 2
  X <- matrix(1, nrow = m, ncol = n)
  w <- Variable(n)
  
  expr2 <- lapply(1:m, function(i) { log_sum_exp(vstack(0, X[i,] %*% w)) })
  expr3 <- Reduce("+", expr2)
  obj <- Minimize(expr3)
  p <- Problem(obj)
  solve(p, solver = "SCS", max_iters = 1)
})
