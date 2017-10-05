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

test_that("Test Huber regression", {
  n <- 1
  m <- 450
  M <- 1      # Huber threshold
  p <- 0.1    # Fraction of responses with sign flipped
  
  # Generate problem data
  beta_true <- 5*matrix(rnorm(n), nrow = n)
  X <- matrix(rnorm(m*n), nrow = m, ncol = n)
  y_true <- X %*% beta_true
  eps <- matrix(rnorm(m), nrow = m)
  
  # Randomly flip sign of some responses
  factor <- 2*rbinom(m, size = 1, prob = 1-p) - 1
  y <- factor * y_true + eps
  
  # Solve ordinary least squares problem
  beta <- Variable(n)
  rel_err <- norm(beta - beta_true, "F")/norm(beta_true, "F")
  
  obj <- sum((y - X %*% beta)^2)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  beta_ols <- result$getValue(beta)
  err_ols <- result$getValue(rel_err)
  
  # Plot fit against measured responses
  plot(X[factor == 1], y[factor == 1], col = "black", xlab = "X", ylab = "y")
  points(X[factor == -1], y[factor == -1], col = "red")
  lines(X, X %*% beta_ols, col = "blue")
  
  # Solve Huber regression problem
  obj <- sum(Huber(y - X %*% beta, M))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  beta_hub <- result$getValue(beta)
  err_hub <- result$getValue(rel_err)
  lines(X, X %*% beta_hub, col = "seagreen", lty = "dashed")
  
  # Solve ordinary least squares assuming sign flips known
  obj <- sum((y - factor*(X %*% beta))^2)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  beta_prs <- result$getValue(beta)
  err_prs <- result$getValue(rel_err)
  lines(X, X %*% beta_prs, col = "black")
  legend("topright", c("OLS", "Huber", "Prescient"), col = c("blue", "seagreen", "black"), lty = 1)
})

test_that("Test logistic regression", {
  n <- 20
  m <- 1000
  offset <- 0
  sigma <- 45
  DENSITY <- 0.2
  
  beta_true <- rnorm(n)
  idxs <- sample(n, size = floor((1-DENSITY)*n), replace = FALSE)
  beta_true[idxs] <- 0
  X <- matrix(rnorm(m*n, 0, 5), nrow = m, ncol = n)
  y <- sign(X %*% beta_true + offset + rnorm(m, 0, sigma))
  
  beta <- Variable(n)
  obj <- sum(Logistic(-X[y <= 0,] %*% beta)) + sum(Logistic(X[y == 1,] %*% beta))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  
  log_odds <- result$getValue(X %*% beta)
  beta_res <- result$getValue(beta)
  y_probs <- 1/(1 + exp(-X %*% beta_res))
  log(y_probs/(1 - y_probs))
})


test_that("Test saturating hinges problem", {
  if(!("ElemStatLearn" %in% rownames(installed.packages())))
    install.packages("ElemStatLearn")
  library(ElemStatLearn)
  
  # Import and sort data
  data(bone)
  X <- bone[bone$gender == "female",]$age
  y <- bone[bone$gender == "female",]$spnbmd
  ord <- order(X, decreasing = FALSE)
  X <- X[ord]
  y <- y[ord]
  
  # Choose knots evenly distributed along domain
  k <- 10
  lambdas <- c(1, 0.5, 0.01)
  idx <- floor(seq(1, length(X), length.out = k))
  knots <- X[idx]
  
  # Saturating hinge
  f_est <- function(x, knots, w0, w) {
    hinges <- sapply(knots, function(t) { pmax(x - t, 0) })
    w0 + hinges %*% w
  }
  
  # Loss function
  loss_obs <- function(y, f) { (y - f)^2 }
  
  # Form problem
  w0 <- Variable(1)
  w <- Variable(k)
  loss <- sum(loss_obs(y, f_est(X, knots, w0, w)))
  constr <- list(sum(w) == 0)
  
  xrange <- seq(min(X), max(X), length.out = 100)
  splines <- matrix(0, nrow = length(xrange), ncol = length(lambdas))
  
  for(i in 1:length(lambdas)) {
    lambda <- lambdas[i]
    reg <- lambda * Norm(w, 1)
    obj <- loss + reg
    prob <- Problem(Minimize(obj), constr)
    
    # Solve problem and save spline weights
    result <- solve(prob)
    w0s <- result$getValue(w0)
    ws <- result$getValue(w)
    splines[,i] <- f_est(xrange, knots, w0s, ws)
  }
  
  # Plot saturating hinges
  plot(X, y, xlab = "Age", ylab = "Change in Bone Density", col = "black", type = "p")
  matlines(xrange, splines, col = "blue", lty = 1:length(lambdas), lwd = 1.5)
  legend("topright", as.expression(lapply(lambdas, function(x) bquote(lambda==.(x)))), col = "blue", lty = 1:length(lambdas), bty = "n")
  
  # Add outliers to data
  set.seed(1)
  nout <- 50
  X_out <- runif(nout, min(X), max(X))
  y_out <- runif(nout, min(y), 3*max(y)) + 0.3
  X_all <- c(X, X_out)
  y_all <- c(y, y_out)
  
  # Solve with squared error loss
  loss_obs <- function(y, f) { (y - f)^2 }
  loss <- sum(loss_obs(y_all, f_est(X_all, knots, w0, w)))
  prob <- Problem(Minimize(loss + reg), constr)
  result <- solve(prob)
  spline_sq <- f_est(xrange, knots, result$getValue(w0), result$getValue(w))
  
  # Solve with Huber loss
  loss_obs <- function(y, f, M) { Huber(y - f, M) }
  loss <- sum(loss_obs(y, f_est(X, knots, w0, w), 0.01))
  prob <- Problem(Minimize(loss + reg), constr)
  result <- solve(prob)
  spline_hub <- f_est(xrange, knots, result$getValue(w0), result$getValue(w))
  
  # Compare fitted functions with squared error and Huber loss
  plot(X, y, xlab = "Age", ylab = "Change in Bone Density", col = "black", type = "p", ylim = c(min(y), 1))
  points(X_out, y_out, col = "red", pch = 16)
  matlines(xrange, cbind(spline_hub, spline_sq), col = "blue", lty = 1:2, lwd = 1.5)
  legend("topright", c("Huber Loss", "Squared Loss"), col = "blue", lty = 1:2, bty = "n")
})

test_that("Test maximum likelihood estimation", {
  set.seed(1)
  
  # Problem data
  m <- 1000
  lambda <- 4
  x <- rpois(m, lambda)
  
  K <- max(x)
  xhist <- hist(x, breaks = -1:K, right = TRUE, include.lowest = FALSE)
  counts <- xhist$counts
  
  # Form problem with log-concave constraint
  u <- Variable(K+1)
  obj <- t(counts) %*% u
  constraints <- list(sum(exp(u)) <= 1, diff(u[1:(K-1)]) >= diff(u[2:K]))
  prob <- Problem(Maximize(obj), constraints)
  result <- solve(prob)
  
  # Plot probability mass function
  pmf <- result$getValue(exp(u))
  plot(0:K, pmf, type = "b", xlab = "x", ylab = "Probability Mass Function")
})

test_that("Test catenary problem", {
  # Problem data
  m <- 101
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
  plot(c(0, 1), c(0, 1), type = 'n', xlab = "x", ylab = "y")
  lines(xs, ys, col = "blue", lwd = 2)
  grid()
  
  ideal <- function(x) { 0.22964*cosh((x-0.5)/0.22964)-0.02603 }
  expect_equal(ys, ideal(xs), tolerance = 1e-3)
  # points(c(0, 1), c(1, 1))
  # curve(0.22964*cosh((x-0.5)/0.22964)-0.02603, 0, 1, col = "red", add = TRUE)
  # grid()
  
  # Lower right endpoint and add staircase structure
  ground <- sapply(seq(0, 1, length.out = m), function(x) {
    if(x < 0.2)
      return(0.6)
    else if(x >= 0.2 && x < 0.4)
      return(0.4)
    else if(x >= 0.4 && x < 0.6)
      return(0.2)
    else
      return(0)
  })
  constraints <- c(constraints, y >= ground)
  constraints[[4]] <- (y[m] == 0.5)
  prob <- Problem(objective, constraints)
  result <- solve(prob)

  # Plot catenary against ground
  xs <- result$getValue(x)
  ys <- result$getValue(y)
  plot(c(0, 1), c(1, 0.5), type = "n", xlab = "x", ylab = "y", ylim = c(0, 1))
  points(c(0, 1), c(1, 0.5))
  lines(xs, ys, col = "blue", lwd = 2)
  lines(xs, ground, col = "red")
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

  # Plot probability density function
  cl <- rainbow(3)
  plot(density(ypop), col = cl[1], xlab = "y", ylab = NA, ylim = c(0, 0.2), zero.line = FALSE)
  lines(density(y), col = cl[2])
  lines(density(y, weights = weights), col = cl[3])
  legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
  
  # Plot cumulative distribution function
  plot(NA, xlab = "y", ylab = NA, xlim = c(-2, 3), ylim = c(0, 1))
  plot_cdf(y, color = cl[1])
  plot_cdf(y[sub], color = cl[2])
  plot_cdf(y[sub], weights, color = cl[3])
  legend("topleft", c("True", "Sample", "Estimate"), lty = c(1,1,1), col = cl)
})

test_that("Test risk-return tradeoff in portfolio optimization", {
  # Problem data
  set.seed(10)
  n <- 10
  SAMPLES <- 100
  mu <- matrix(abs(rnorm(n)), nrow = n)
  Sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
  Sigma <- t(Sigma) %*% Sigma
  
  # Form problem
  w <- Variable(n)
  ret <- t(mu) %*% w
  risk <- QuadForm(w, Sigma)
  constraints <- list(w >= 0, sum(w) == 1)
  
  # Risk aversion parameters
  gammas <- 10^seq(-2, 3, length.out = SAMPLES)
  ret_data <- rep(0, SAMPLES)
  risk_data <- rep(0, SAMPLES)
  w_data <- matrix(0, nrow = SAMPLES, ncol = n)
  
  # Compute trade-off curve
  for(i in 1:length(gammas)) {
    gamma <- gammas[i]
    objective <- ret - gamma * risk
    prob <- Problem(Maximize(objective), constraints)
    result <- solve(prob)
    
    # Evaluate risk/return for current solution
    risk_data[i] <- result$getValue(sqrt(risk))
    ret_data[i] <- result$getValue(ret)
    w_data[i,] <- result$getValue(w)
  }
   
  # Plot trade-off curve
  plot(risk_data, ret_data, xlab = "Risk (Standard Deviation)", ylab = "Return", xlim = c(0, 4), ylim = c(0, 2), type = "l", lwd = 2, col = "blue")
  points(sqrt(diag(Sigma)), mu, col = "red", cex = 1.5, pch = 16)
  markers_on <- c(10, 20, 30, 40)
  for(marker in markers_on) {
    points(risk_data[marker], ret_data[marker], col = "black", cex = 1.5, pch = 15)
    nstr <- sprintf("%.2f", gammas[marker])
    text(risk_data[marker] + 0.2, ret_data[marker] - 0.05, bquote(paste(gamma, " = ", .(nstr))), cex = 1.5)
  }
  
  # Plot weights for a few gamma
  w_plot <- t(w_data[markers_on,])
  colnames(w_plot) <- sprintf("%.2f", gammas[markers_on])
  if("colorspace" %in% rownames(installed.packages())) {
    require(colorspace)
    barplot(w_plot, xlab = expression(paste("Risk Aversion (", gamma, ")", sep = "")), ylab = "Fraction of Budget", col = sequential_hcl(n))
  } else
    barplot(w_plot, xlab = expression(paste("Risk Aversion (", gamma, ")", sep = "")), ylab = "Fraction of Budget")
})

test_that("Test Kelly gambling optimal bets", {
  set.seed(1)
  n <- 20      # Total bets
  K <- 100     # Number of possible returns
  PERIODS <- 100
  TRIALS <- 5
  
  # Generate return probabilities
  ps <- runif(K)
  ps <- ps/sum(ps)
  
  # Generate matrix of possible returns
  rets <- runif(K*(n-1), 0.5, 1.5)
  shuff <- sample(1:length(rets), size = length(rets), replace = FALSE)
  rets[shuff[1:30]] <- 0    # Set 30 returns to be relatively low
  rets[shuff[31:60]] <- 5     # Set 30 returns to be relatively high
  rets <- matrix(rets, nrow = K, ncol = n-1)
  rets <- cbind(rets, rep(1, K))   # Last column represents not betting
    
  # Solve for Kelly optimal bets
  b <- Variable(n)
  obj <- Maximize(t(ps) %*% log(rets %*% b))
  constraints <- list(sum(b) == 1, b >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  bets <- result$getValue(b)

  # Naive betting scheme: equal split on bets with highest expected returns
  bets_cmp <- matrix(0, nrow = n)
  bets_cmp[n] <- 0.15                  # Hold 15% of wealth
  rets_avg <- ps %*% rets
  tidx <- order(rets_avg[-n], decreasing = TRUE)[1:9]
  fracs <- rets_avg[tidx]/sum(rets_avg[tidx])
  bets_cmp[tidx] <- fracs*(1-bets_cmp[n])
  
  # Calculate wealth over time
  wealth <- matrix(0, nrow = PERIODS, ncol = TRIALS)
  wealth_cmp <- matrix(0, nrow = PERIODS, ncol = TRIALS)
  for(i in 1:TRIALS) {
    sidx <- sample(1:K, size = PERIODS, replace = TRUE, prob = ps)
    winnings <- rets[sidx,] %*% bets
    wealth[,i] <- cumprod(winnings)
    
    winnings_cmp <- rets[sidx,] %*% bets_cmp
    wealth_cmp[,i] <- cumprod(winnings_cmp)
  }
  
  # Plot Kelly optimal growth trajectories
  matplot(1:PERIODS, wealth, xlab = "Time", ylab = "Wealth", log = "y", type = "l", col = "blue", lty = 1, lwd = 2)
  matlines(1:PERIODS, wealth_cmp, col = "red", lty = 2, lwd = 2)
  legend("topleft", c("Kelly Optimal Bets", "Naive Bets"), col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
})

test_that("Test worst-case covariance", {
  # Problem data
  w <- matrix(c(0.1, 0.2, -0.05, 0.1))
  # Constraint matrix:
  # [[0.2, + ,  +,   +-],
  #  [+,   0.1, -,   - ],
  #  [+,   -,   0.3, + ],
  #  [+-,  -,   +,  0.1]]

  # Form problem
  Sigma <- Semidef(4)
  obj <- Maximize(t(w) %*% Sigma %*% w)
  constraints <- list(Sigma[1,1] == 0.2, Sigma[2,2] == 0.1, Sigma[3,3] == 0.3, Sigma[4,4] == 0.1,
                      Sigma[1,2] >= 0, Sigma[1,3] >= 0, Sigma[2,3] <= 0, Sigma[2,4] <= 0, Sigma[3,4] >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  
  Sigma_true <- rbind(c(0.2,    0.0978, 0,     0.0741),
                      c(0.0978, 0.1,   -0.101, 0),
                      c(0,     -0.101,  0.3,   0),
                      c(0.0741, 0,      0,     0.1))
  dimnames(Sigma_true) <- list(NULL, NULL)
  expect_equal(sqrt(result$value), 0.1232, tolerance = 1e-3)
  expect_equal(as.matrix(result$getValue(Sigma)), Sigma_true, tolerance = 1e-3)
  
  # Problem data
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
  constr <- list(L <= Sigma, Sigma <= U, Sigma == t(Sigma))
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob)
  print(result$getValue(Sigma))
})

test_that("Test sparse inverse covariance estimation", {
  require(Matrix)
  require(expm)
  
  set.seed(1)
  n <- 10      # Dimension of matrix
  m <- 1000    # Number of samples
  alphas <- c(10, 4, 1)
  
  # Create sparse, symmetric PSD matrix S
  A <- rsparsematrix(n, n, 0.15, rand.x = rnorm)
  Strue <- A %*% t(A) + 0.05 * diag(rep(1, n))    # Force matrix to be strictly positive definite
  image(Strue != 0, main = "Inverse of true covariance matrix")
  
  # Create covariance matrix associated with S
  R <- base::solve(Strue)
  
  # Sample from distribution with covariance R
  # If Y ~ N(0, I), then R^{1/2} * Y ~ N(0, R) since R is symmetric
  x_sample <- matrix(rnorm(n*m), nrow = m, ncol = n) %*% t(expm::sqrtm(R))
  Q <- cov(x_sample)    # Sample covariance matrix
  
  S <- Semidef(n)    # Variable constrained to positive semidefinite cone
  obj <- Maximize(LogDet(S) - Trace(S %*% Q))
  for(alpha in alphas) {
    constraints <- list(sum(abs(S)) <= alpha)
  
    # Form and solve optimization problem
    prob <- Problem(obj, constraints)
    result <- solve(prob)
  
    # Create covariance matrix
    R_hat <- base::solve(result$getValue(S))
    Sres <- result$getValue(S)
    Sres[abs(Sres) <= 1e-4] <- 0
    title <- bquote(bold(paste("Estimated inv. cov matrix (", alpha, " = ", .(alpha), ")")))
    image(Sres != 0, main = title)
  }
})

test_that("Test fastest mixing Markov chain (FMMC)", {
  # Boyd, Diaconis, and Xiao. SIAM Rev. 46 (2004) pgs. 667-689 at pg. 672
  # Form the complementary graph
  antiadjacency <- function(g) {
    n <- max(as.numeric(names(g)))   # Assumes names are integers starting from 1
    a <- lapply(1:n, function(i) c())
    names(a) <- 1:n
    for(x in names(g)) {
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
    for(i in names(a)) {
      for(j in a[[i]]) {  # (i-j) is a not-edge of g!
        idx <- as.numeric(i)
        if(idx != j)
          constraints <- c(constraints, P[idx,j] == 0)
      }
    }
    prob <- Problem(objective, constraints)
    result <- solve(prob)
    if(verbose)
      cat("Status: ", result$status, ", Optimal Value = ", result$value)
    list(status = result$status, value = result$value, P = result$getValue(P))
  }
  
  disp_result <- function(states, P, tol = 1e-3) {
    if(!("markovchain" %in% rownames(installed.packages()))) {
      rownames(P) <- states
      colnames(P) <- states
      print(P)
    } else {
      require(markovchain)
      P[P < tol] <- 0
      P <- P/apply(P, 1, sum)   # Normalize so rows sum to exactly 1
      mc <- new("markovchain", states = states, transitionMatrix = P)
      plot(mc, label.cex = 1.5)
    }
  }
  
  # SIAM Rev. 46 examples pg. 674: Figure 1 and Table 1
  # a) line graph L(4)
  g <- list("1" = 2, "2" = c(1,3), "3" = c(2,4), "4" = 3)
  result <- FMMC(g, verbose = TRUE)
  disp_result(names(g), result$P)
  
  # b) triangle + one edge
  g <- list("1" = 2, "2" = c(1,3,4), "3" = c(2,4), "4" = c(2,3))
  result <- FMMC(g, verbose = TRUE)
  disp_result(names(g), result$P)
  
  # c) bipartite 2 + 3
  g <- list("1" = c(2,4,5), "2" = c(1,3), "3" = c(2,4,5), "4" = c(1,3), "5" = c(1,3))
  result <- FMMC(g, verbose = TRUE)
  disp_result(names(g), result$P)
  
  # d) square + central point
  g <- list("1" = c(2,3,5), "2" = c(1,4,5), "3" = c(1,4,5), "4" = c(2,3,5), "5" = c(1,2,3,4,5))
  result <- FMMC(g, verbose = TRUE)
  disp_result(names(g), result$P)
})
