test_that("Test non-negative least squares", {
  require(MASS)
  
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
  objective <- Minimize(sum((y - X %*% beta)^2))
  prob <- Problem(objective)
  
  # Solve the OLS problem for beta
  system.time(result <- solve(prob))
  beta_ols <- result$getValue(beta)
  
  # Add non-negativity constraint on beta
  constraints <- list(beta >= 0)
  prob2 <- Problem(objective, constraints)
  
  # Solve the NNLS problem for beta
  system.time(result2 <- solve(prob2))
  beta_nnls <- result$getValue(beta)
  expect_true(all(beta_nnls >= 0))   # All resulting beta should be non-negative
  
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
  m <- 51
  L <- 2
  h <- L/(m-1)
  
  # Form objective
  x <- Variable(m)
  y <- Variable(m)
  objective <- Minimize(sum(y))
  
  # Form constraints
  constraints <- list(x[1] == 0, y[1] == 1, x[m] == 1, y[m] == 1, 
                      diff(x)^2 + diff(y)^2 <= h^2)
  
  # Solve the catenary problem
  prob <- Problem(objective, constraints)
  system.time(result <- solve(prob))
  
  # Plot and compare with ideal catenary
  xs <- result$getValue(x)
  ys <- result$getValue(y)
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
  objective <- sum(Entr(w))
  constraints <- list(w >= 0, sum(w) == 1, t(X[sub,]) %*% w == b)
  prob <- Problem(Maximize(objective), constraints)
  
  # Solve for the distribution weights
  result <- solve(prob)
  weights <- result$getValue(w)

  cl <- rainbow(3)
  plot(NA, xlab = "y", ylab = NA, xlim = c(-2, 3), ylim = c(0, 1))
  plot_cdf(y, color = cl[1])
  plot_cdf(y[sub], color = cl[2])
  plot_cdf(y[sub], weights, color = cl[3])
  legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
})

test_that("Test risk-return trade-off in portfolio optimization", {
  # Problem data
  n <- 10
  mu <- matrix(abs(rnorm(n)), nrow = n)
  Sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
  Sigma <- t(Sigma) %*% Sigma
  
  # Form problem
  w <- Variable(n)
  ret <- t(mu) %*% w
  risk <- QuadForm(w, Sigma)
  constraints <- list(w >= 0, sum(w) == 1)
  
  # Risk aversion parameters
  SAMPLES <- 100
  gammas <- 10^seq(-2, 3, length.out = SAMPLES)
  ret_data <- rep(0, SAMPLES)
  risk_data <- rep(0, SAMPLES)
  
  # Compute trade-off curve
  for(gamma in gammas) {
    objective <- ret - gamma * risk
    prob <- Problem(Maximize(objective), constraints)
    result <- solve(prob)
    
    # Evaluate risk/return for current solution
    ret_data[i] <- result$getValue(ret)
    risk_data[i] <- result$getValue(risk)
  }
   
  # Plot trade-off curve
  plot(risk, return, main = "Risk-Return Curve", xlab = "Variance", ylab = "Return")
  points(diag(Sigma), mu, col = "red")
})

test_that("Test worst-case covariance", {
  n <- 5
  w <- rexp(n)
  w <- w / sum(w)
  mu <- abs(matrix(rnorm(n), nrow = n, ncol = 1))/15
  Sigma_nom <- matrix(runif(n^2, -0.15, 0.8), nrow = n, ncol = n)
  Sigma_nom <- t(Sigma_nom) %*% Sigma_nom
  
  # Known upper and lower bounds
  Delta <- matrix(0.2, nrow = n, ncol = n)
  diag(Delta) <- 0
  U <- Sigma_nom + Delta
  L <- Sigma_nom - Delta
  
  Sigma <- Semidef(n)
  obj <- QuadForm(w, Sigma)
  constr <- list(L <= Sigma, Sigma <= U, Sigma == Symmetric(n))
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob)
  print(result$getValue(Sigma))
})

test_that("Test fastest mixing Markov chain (FMMC)", {
  # Boyd, Diaconis, and Xiao. SIAM Rev. 46 (2004) pgs. 667-689 at pg. 672
  
  # Form the complementary graph
  antiadjacency <- function(g) {
    n <- max(as.numeric(names(g)))   # Assumes names are integers starting from 1
    a <- lapply(1:n, function(i) c())
    names(a) <- 1:n
    for(x in g) {
      for(y in 1:n) {
        if(!(y %in% g[[x]]))
          a[[x]] <- c(a[[x]], y)
      }
    }
    a
  }
  
  # Fastest mixing Markov chain on graph g
  FMMC <- function(g, verbose = FALSE) {
    a <- antiadjacency(g)
    n <- length(names(a))
    P <- Variable(n, n)
    o <- rep(1, n)
    objective <- Minimize(Norm(P - 1.0/n))
    constraints <- list(P %*% o == o, t(P) == P, P >= 0)
    for(i in a) {
      for(j in a[[i]]) {  # (i-j) is a not-edge of g!
        if(i != j)
          constraints <- c(constraints, P[i,j] == 0)
      }
    }
    prob <- Problem(objective, constraints)
    result <- solve(prob)
    if(verbose)
      cat("Status: ", result$status, ", Optimal Value = ", result$value)
    list(status = result$status, value = result$value, P = result$getValue(P))
  }
  
  # SIAM Rev. 46 examples pg. 674: Figure 1 and Table 1
  # a) line graph L(4)
  g <- list("1" = 2, "2" = c(1,3), "3" = c(2,4), "4" = 3)
  result <- FMMC(g, verbose = TRUE)
  print(result$P)
  
  # b) triangle + one edge
  g <- list("1" = 2, "2" = c(1,3,4), "3" = c(2,4), "4" = c(2,3))
  result <- FMMC(g, verbose = TRUE)
  print(result$P)
  
  # c) bipartite 2 + 3
  g <- list("1" = c(2,4,5), "2" = c(1,3), "3" = c(2,4,5), "4" = c(1,3), "5" = c(1,3))
  result <- FMMC(g, verbose = TRUE)
  print(result$P)
  
  # d) square + central point
  g <- list("1" = c(2,3,5), "1" = c(1,4,5), "2" = c(1,4,5), "3" = c(2,3,5), "4" = c(1,2,3,4,5))
  result <- FMMC(g, verbose = TRUE)
  print(result$P)
})
