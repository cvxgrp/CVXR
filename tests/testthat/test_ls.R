TOL <- 1e-6

test_that("Test regression", {
  # Set the random seed to get consistent data
  set.seed(1)
  
  # Number of examples to use
  n <- 100
  
  # Specify the true value of the variable
  true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)
  
  # Generate data
  x_data <- matrix(runif(n) * 5, nrow = n, ncol = 1)
  x_data_expanded <- cbind(x_data, x_data^2, x_data^3)
  y_data <- x_data_expanded %*% true_coeffs + 0.5 * matrix(runif(n, 1), nrow = n, ncol = 1)
  
  slope <- Variable()
  offset <- Variable()
  line <- offset + x_data * slope
  residuals <- line - y_data
  fit_error <- SumSquares(residuals)
  # result <- solve(Problem(Minimize(fit_error), list()), solver = "LS")
  # optval <- result$optimal_value
  # expect_equal(optval, 1171.60037715, tolerance = TOL)
  
  quadratic_coeff <- Variable()
  slope <- Variable()
  offset <- Variable()
  quadratic <- offset + x_data %*% slope + quadratic_coeff %*% x_data^2
  residuals <- quadratic - y_data
  fit_error <- SumSquares(residuals)
  # result <- solve(Problem(Minimize(fit_error), list()), solver = "LS")
  result2 <- solve(Problem(Minimize(fit_error), list()), solver = "ECOS")
  # optval <- result$optimal_value
  # optval2 <- result2$optimal_value
  # expect_equal(optval, 139.225650756, tolerance = TOL)
})

test_that("Test control", {
  # Some constraints on our motion
  # The object should start from the origin and end at rest
  initial_velocity <- matrix(c(-20, 100), nrow = 2, ncol = 1)
  final_position <- matrix(c(100, 100), nrow = 2, ncol = 1)
  
  T <- 100   # The number of timesteps
  h <- 0.1   # The time between time intervals
  mass <- 1  # Mass of object
  drag <- 0.1  # Drag on object
  g <- matrix(c(0, -9.8), nrow = 2, ncol = 1)   # Gravity on object
  
  # Declare the variables we need
  position <- Variable(2, T)
  velocity <- Variable(2, T)
  force <- Variable(2, T-1)
  
  # Create a problem instance
  mu <- 1
  constraints <- list()
  
  # Add constraints on our variables
  for(i in 1:(T-1)) {
    constraints <- c(constraints, position[,i+1] == position[,i] + h * velocity[,i])
    acceleration <- force[,i]/mass + g - drag * velocity[,i]
    constraints <- c(constraints, velocity[,i+1] == velocity[,i] + h * acceleration)
  }
  
  # Add position constraints
  constraints <- c(constraints, position[,1] == 0)
  constraints <- c(constraints, position[,T] == final_position)
  
  # Add velocity constraints
  constraints <- c(constraints, velocity[,1] == initial_velocity)
  constraints <- c(constraints, velocity[,T] == 0)
  
  # Solve the problem
  # result <- solve(Problem(Minimize(SumSquares(force)), constraints), solver = "LS")
  # optval <- result$optimal_value
  # expect_equal(optval, 17859.0, tolerance = 1)
})

test_that("Test sparse system", {
  m <- 1000
  n <- 800
  r <- 700
  
  set.seed(1)
  density <- 0.2
  A <- rsparsematrix(m, n, density, rand.x = runif)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  G <- rsparsematrix(r, n, density, rand.x = runif)
  h <- matrix(rnorm(r), nrow = r, ncol = 1)
  
  x <- Variable(n)
  # result <- solve(Problem(Minimize(SumSquares(A*x - b)), list(G*x == h)), solver = "LS")
  # optval <- result$optimal_value
  # expect_equal(optval, 6071.830658, tolerance = TOL)
})

test_that("Test equivalent forms", {
  m <- 100
  n <- 80
  r <- 70
  
  set.seed(1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  G <- matrix(rnorm(r*n), nrow = r, ncol = n)
  h <- matrix(rnorm(r), nrow = r, ncol = 1)
  
  # ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2
  P <- t(A) %*% A
  q <- -2 * t(A) %*% b
  r <- t(b) %*% b
  
  Pinv <- solve(P)
  
  x <- Variable(n)
  
  obj1 <- SumSquares(A %*% x - b)
  obj2 <- SumEntries(Square(A %*% x - b))
  obj3 <- QuadForm(x, P) + t(q) %*% x + r
  obj4 <- MatrixFrac(x, Pinv) + t(x) %*% x + r
  
  cons <- list(G %*% x == h)
  
  # v1 <- solve(Problem(Minimize(obj1), cons), solver = "LS")$optimal_value
  # v2 <- solve(Problem(Minimize(obj2), cons), solver = "LS")$optimal_value
  # v3 <- solve(Problem(Minimize(obj3), cons), solver = "LS")$optimal_value
  # v4 <- solve(Problem(Minimize(obj4), cons), solver = "LS")$optimal_value
  
  # expect_equal(v1, 681.119420108, tolerance = TOL)
  # expect_equal(v2, 681.119420108, tolerance = TOL)
  # expect_equal(v3, 681.119420108, tolerance = TOL)
  # expect_equal(v4, 681.119420108, tolerance = TOL)
})

test_that("Test smooth ridge", {
  set.seed(1)
  n <- 500
  k <- 50
  delta <- 1
  eta <- 1
  
  A <- matrix(runif(k*n), nrow = k, ncol = n)
  b <- matrix(runif(k), nrow = k, ncol = 1)
  x <- Variable(n)
  obj <- SumSquares(A %*% x - b) + delta*SumSquares(x[1:(n-1)]-x[2:n]) + eta*SumSquares(x)
  # optval <- solve(Problem(Minimize(obj), list()), solver = "LS")$optimal_value
  # expect_equal(optval, 0.24989717371, tolerance = TOL)
})
