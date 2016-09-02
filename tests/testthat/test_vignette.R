test_that("Test non-negative least squares", {
  library(MASS)
  
  # Generate problem data
  s <- 1
  m <- 10
  n <- 300
  mu <- rep(0, 9)
  Sigma <- data.frame(c(1.6484, -0.2096, -0.0771, -0.4088, 0.0678, -0.6337, 0.9720, -1.2158, -1.3219),
                      c(-0.2096, 1.9274, 0.7059, 1.3051, 0.4479, 0.7384, -0.6342, 1.4291, -0.4723),
                      c(-0.0771, 0.7059, 2.5503, 0.9047, 0.9280, 0.0566, -2.5292, 0.4776, -0.4552),
                      c(-0.4088, 1.3051, 0.9047, 2.7638, 0.7607, 1.2465, -1.8116, 2.0076, -0.3377),
                      c(0.0678, 0.4479, 0.9280, 0.7607, 3.8453, -0.2098, -2.0078, -0.1715, -0.3952),
                      c(-0.6337, 0.7384, 0.0566, 1.2465, -0.2098, 2.0432, -1.0666,  1.7536, -0.1845),
                      c(0.9720, -0.6342, -2.5292, -1.8116, -2.0078, -1.0666, 4.0882,  -1.3587, 0.7287),
                      c(-1.2158, 1.4291, 0.4776, 2.0076, -0.1715, 1.7536, -1.3587, 2.8789, 0.4094),
                      c(-1.3219, -0.4723, -0.4552, -0.3377, -0.3952, -0.1845, 0.7287, 0.4094, 4.8406))
  
  X <- mvrnorm(n, mu, Sigma)
  X <- cbind(rep(1, n), X)
  b <- c(0, 0.8, 0, 1, 0.2, 0, 0.4, 1, 0, 0.7)
  y <- X %*% b + rnorm(n, 0, s)
  
  # Construct the OLS problem without constraints
  beta <- Variable(m)
  objective <- Minimize(SumSquares(y - X %*% beta))
  prob <- Problem(objective)
  
  # Solve the OLS problem for beta
  system.time(result <- cvxr_solve(prob))
  result$optimal_value
  result$primal_values[[as.character(beta@id)]]
  beta_ols <- result$primal_values[[as.character(beta@id)]]
  
  # Add non-negativity constraint on beta
  constraints <- list(beta >= 0)
  prob2 <- Problem(objective, constraints)
  
  # Solve the NNLS problem for beta
  system.time(result2 <- cvxr_solve(prob2))
  result2$optimal_value
  result2$primal_values[[as.character(beta@id)]]
  beta_nnls <- result2$primal_values[[as.character(beta@id)]]
  all(beta_nnls >= 0)   # All resulting beta should be non-negative
  
  # Calculate the fitted y values
  fit_ols <- X %*% beta_ols
  fit_nnls <- X %*% beta_nnls
  
  # Plot coefficients for OLS and NNLS
  coeff <- cbind(b, beta_ols, beta_nnls)
  colnames(coeff) <- c("Actual", "OLS", "NNLS")
  rownames(coeff) <- paste("beta", 1:length(b)-1, sep = "")
  barplot(t(coeff), ylab = "Coefficients", beside = TRUE, legend = TRUE)
})

test_that("Test catenary problem", {
  # Problem data
  n <- 51
  L <- 2
  h <- L/(n-1)
  
  # Form objective
  x <- Variable(2*n)
  B <- diag(2*n)
  B[1:n, 1:n] <- 0
  objective <- Minimize(SumEntries(B %*% x))
  
  # Form constraints
  A <- matrix(0, nrow = 4, ncol = 2*n)
  A[1, 1] <- A[2, n] <- A[3, n+1] <- A[4, 2*n] <- 1
  b <- matrix(c(0, 1, 1, 1), nrow=4)
  constraints = list(x >= 0, A %*% x == b)
  
  for (i in seq.int(n-1)) {
    A <- matrix(numeric(2*2*n), nrow = 2)
    A[1, i] <- -1; A[1, i+1] <- 1
    A[2, n+i] <- -1; A[2, n+i+1] <- 1
    constraints <- c(constraints, Norm2(A %*% x) <= h)
  }
  
  # Solve the catenary problem
  prob <- Problem(objective, constraints)
  system.time(result <- cvxr_solve(prob))
  
  # Plot results
  x <- result$primal_values[[as.character(x@id)]]
  xs <- x[1:n, 1, drop = TRUE]
  ys <- x[(n+1):(2*n), 1, drop = TRUE]
  plot(c(0, 1), c(0, 1), type = 'n')
  lines(xs, ys, col = "blue", lwd = 2)
  
  points(c(0, 1), c(1, 1))
  curve(0.22964*cosh((x-0.5)/0.22964)-0.02603, 0, 1, col = "red", add = TRUE)
  grid()
})

test_that("Test direct standardization problem", {
  skew_sample <- function(data, bias) {
    if(missing(bias))
      bias <- rep(1.0, ncol(data))
    num <- exp(data %*% bias)
    return(num / sum(num))
  }
  
  plot_cdf <- function(data, probs, color = 'k') {
    if(missing(probs))
      probs <- rep(1.0/length(data), length(data))
    distro <- cbind(data, probs)
    dsort <- distro[order(distro[,1]),]
    ecdf <- cumsum(dsort[,2])
    lines(dsort[,1], ecdf, col = color)
  }
  
  # Problem data
  n <- 2
  m <- 1000
  msub <- 100
  
  # Generate original distribution
  sex <- rbinom(m, 1, 0.5)
  age <- sample(10:60, m, replace = TRUE)
  mu <- 5 * sex + 0.1 * age
  X <- cbind(sex, age)
  y <- rnorm(mu, 1.0)
  b <- as.matrix(apply(X, 2, mean))
  
  # Generate skewed subsample
  skew <- skew_sample(X, c(-0.95, -0.05))
  sub <- sample(1:m, msub, replace = TRUE, prob = skew)
  
  # Construct the direct standardization problem
  w <- Variable(msub)
  objective <- SumEntries(Entr(w))
  constraints <- list(w >= 0, SumEntries(w) == 1, t(X[sub,]) %*% w == b)
  prob <- Problem(Maximize(objective), constraints)
  
  # Solve for the distribution weights
  result <- cvxr_solve(prob)
  result$optimal_value
  result$primal_values[[as.character(w@id)]]
  weights <- result$primal_values[[as.character(w@id)]]
  
  cl <- rainbow(3)
  plot(NA, xlab = "y", ylab = NA, xlim = c(-2, 3), ylim = c(0, 1))
  plot_cdf(y, color = cl[1])
  plot_cdf(y[sub], color = cl[2])
  plot_cdf(y[sub], weights, color = cl[3])
  legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
})

test_that("Test risk-return trade-off in portfolio optimization", {
  pbar <- matrix(c(0.12, 0.10, 0.07, 0.03), nrow = 4, ncol = 1)
  sigma <- rbind(c( 0.0064, 0.0008, -0.0011, 0),
                 c( 0.0008, 0.0025,  0,      0),
                 c(-0.0011, 0,       0.0004, 0),
                 c( 0,      0,       0,      0))
  
  x <- Variable(4)
  mu <- Parameter(sign = "positive")
  objective <- -t(pbar) %*% x + mu * t(x) %*% sigma %*% x
  constraints <- list(SumEntries(x) == 1, x >= 0)
  prob <- Problem(Minimize(objective), constraints)
   
  # result <- cvxr_solve(prob)
  # result$optimal_value
  # result$primal_values
})
