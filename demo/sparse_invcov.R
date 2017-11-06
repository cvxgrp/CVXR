if(!require(Matrix))
  stop("Please install the Matrix library")
if(!require(expm))
  stop("Please install the expm library")

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
obj <- Maximize(log_det(S) - matrix_trace(S %*% Q))
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
