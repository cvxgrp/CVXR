## Integration tests derived from cvxr_docs examples
## Each test is self-contained with its own data generation and expected values.

# ═══════════════════════════════════════════════════════════════════
# 1. intro.Rmd — Least squares (LP/QP)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("intro: OLS via CVXR matches lm()", {
  skip_if_not_installed("clarabel")
  set.seed(123)
  n <- 100; p <- 10
  beta <- -4:5
  X <- matrix(rnorm(n * p), nrow = n)
  Y <- X %*% beta + rnorm(n)
  betaHat <- Variable(p)
  objective <- Minimize(sum((Y - X %*% betaHat)^2))
  problem <- Problem(objective)
  psolve(problem)
  expect_equal(as.numeric(value(betaHat)), unname(coef(lm(Y ~ 0 + X))),
               tolerance = 1e-4)
})

## @cvxpy NONE
test_that("intro: nonneg least squares", {
  skip_if_not_installed("clarabel")
  set.seed(123)
  n <- 100; p <- 10
  beta <- -4:5
  X <- matrix(rnorm(n * p), nrow = n)
  Y <- X %*% beta + rnorm(n)
  betaHat <- Variable(p)
  objective <- Minimize(sum((Y - X %*% betaHat)^2))
  problem <- Problem(objective, constraints = list(betaHat >= 0))
  psolve(problem)
  sol <- as.numeric(value(betaHat))
  expect_true(all(sol >= -1e-6))
})

# ═══════════════════════════════════════════════════════════════════
# 2. portfolio-optimization.Rmd — QP (Markowitz)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("portfolio: Markowitz mean-variance", {
  skip_if_not_installed("clarabel")
  set.seed(10)
  n <- 10
  mu <- abs(rnorm(n, mean = 0.03))
  Sigma_factor <- matrix(rnorm(n^2), nrow = n)
  Sigma <- t(Sigma_factor) %*% Sigma_factor
  w <- Variable(n)
  ret <- t(mu) %*% w
  risk <- quad_form(w, Sigma)
  gamma_val <- 1.0
  prob <- Problem(Maximize(ret - gamma_val * risk),
                  list(w >= 0, sum(w) == 1))
  psolve(prob)
  w_val <- as.numeric(value(w))
  expect_true(all(w_val >= -1e-6))
  expect_equal(sum(w_val), 1, tolerance = 1e-6)
  expect_true(as.numeric(value(ret)) > 0)
})

# ═══════════════════════════════════════════════════════════════════
# 3. huber-regression.Rmd — SOCP (Huber atom)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("huber regression: Huber beats OLS on corrupted data", {
  skip_if_not_installed("clarabel")
  set.seed(1289)
  n <- 1; m <- 450
  beta_true <- 5 * rnorm(n)
  X <- matrix(rnorm(m * n), nrow = m)
  y_true <- X %*% beta_true
  eps <- rnorm(m, 0, 1)
  y <- y_true + eps
  ## Corrupt 10% of observations
  p_corrupt <- 0.1
  idx <- sample.int(m, size = ceiling(m * p_corrupt))
  y[idx] <- -y[idx]

  beta <- Variable(n)
  ## OLS
  prob_ols <- Problem(Minimize(sum((y - X %*% beta)^2)))
  psolve(prob_ols)
  beta_ols <- as.numeric(value(beta))

  ## Huber
  M <- 1
  prob_hub <- Problem(Minimize(sum(huber(y - X %*% beta, M))))
  psolve(prob_hub)
  beta_hub <- as.numeric(value(beta))

  ## Huber should be closer to the true beta
  err_ols <- abs(beta_ols - beta_true) / abs(beta_true)
  err_hub <- abs(beta_hub - beta_true) / abs(beta_true)
  expect_true(err_hub < err_ols)
})

# ═══════════════════════════════════════════════════════════════════
# 4. catenary.Rmd — SOCP (hanging chain)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("catenary: hanging chain problem", {
  skip_if_not_installed("clarabel")
  m <- 101
  L <- 2; h <- L / (m - 1)
  x <- Variable(m)
  y <- Variable(m)
  objective <- Minimize(sum(y))
  constraints <- list(
    x[1] == 0, y[1] == 1,
    x[m] == 1, y[m] == 1,
    diff(x)^2 + diff(y)^2 <= h^2
  )
  prob <- Problem(objective, constraints)
  opt_val <- psolve(prob)
  ys <- as.numeric(value(y))
  ## Chain should sag below endpoints
  expect_true(min(ys) < 1)
  ## Endpoints fixed
  expect_equal(ys[1], 1, tolerance = 1e-4)
  expect_equal(ys[m], 1, tolerance = 1e-4)
  ## Objective should be finite and reasonable
  expect_true(is.finite(opt_val))
  expect_true(opt_val < m)  # less than if all y=1
})

# ═══════════════════════════════════════════════════════════════════
# 5. logistic-regression.Rmd — ExpCone (logistic atom)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("logistic regression: logistic atom fits correctly", {
  skip_if_not_installed("clarabel")
  set.seed(183991)
  n <- 5; m <- 100  # smaller than docs for speed
  beta_true <- matrix(rnorm(n), nrow = n)
  X <- matrix(rnorm(m * n, 0, 5), nrow = m)
  y <- sign(X %*% beta_true + rnorm(m, 0, 5))

  beta <- Variable(n)
  ## y ∈ {0, 1} formulation
  y01 <- (y + 1) / 2
  obj <- Maximize(-sum(X[y <= 0, ] %*% beta) - sum(logistic(-X %*% beta)))
  prob <- Problem(obj)
  psolve(prob)
  beta_val <- as.numeric(value(beta))
  ## Predictions should correlate with true labels
  preds <- sign(X %*% beta_val)
  accuracy <- mean(preds == y)
  expect_true(accuracy > 0.6)
})

# ═══════════════════════════════════════════════════════════════════
# 6. elastic-net.Rmd — QP + L1 (elastic net)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("elastic net: regularization shrinks coefficients", {
  skip_if_not_installed("clarabel")
  set.seed(1)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p, sd = 5), nrow = n)
  beta_true <- c(rep(5, 3), rep(0, 7))
  y <- X %*% beta_true + rnorm(n, sd = 1)

  beta <- Variable(p)
  loss <- sum((y - X %*% beta)^2) / (2 * n)
  alpha <- 0.75
  lambda_small <- 0.01
  lambda_large <- 10

  ## Small lambda: close to OLS
  elastic_small <- lambda_small * ((1 - alpha) / 2 * sum(beta^2) +
                                     alpha * p_norm(beta, 1))
  prob_small <- Problem(Minimize(loss + elastic_small))
  psolve(prob_small)
  beta_small <- as.numeric(value(beta))

  ## Large lambda: heavily shrunk
  elastic_large <- lambda_large * ((1 - alpha) / 2 * sum(beta^2) +
                                     alpha * p_norm(beta, 1))
  prob_large <- Problem(Minimize(loss + elastic_large))
  psolve(prob_large)
  beta_large <- as.numeric(value(beta))

  ## Large lambda should have smaller coefficients
  expect_true(sum(abs(beta_large)) < sum(abs(beta_small)))
})

# ═══════════════════════════════════════════════════════════════════
# 7. kelly-strategy.Rmd — ExpCone (log utility)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("kelly: log-optimal betting", {
  skip_if_not_installed("clarabel")
  set.seed(1)
  n_bets <- 5
  K <- 20  # scenarios (smaller than docs for speed)
  ps <- runif(K)
  ps <- ps / sum(ps)
  rets <- matrix(runif(K * (n_bets - 1), 0.5, 1.5), nrow = K)
  rets <- cbind(rets, rep(1, K))  # no-bet option
  n <- ncol(rets)

  b <- Variable(n)
  obj <- Maximize(t(ps) %*% log(rets %*% b))
  prob <- Problem(obj, list(sum(b) == 1, b >= 0))
  psolve(prob)
  bets <- as.numeric(value(b))
  expect_true(all(bets >= -1e-6))
  expect_equal(sum(bets), 1, tolerance = 1e-6)
})

# ═══════════════════════════════════════════════════════════════════
# 8. integer-programming.Rmd — MILP
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("integer programming: simple MILP", {
  skip_if_not_installed("highs")
  y1 <- Variable(2)
  y2 <- Variable(1, integer = TRUE)
  y3 <- Variable(1)
  x <- vstack(y1, y2, y3)
  C <- matrix(c(1, 2, -0.1, -3), nrow = 1)
  objective <- Maximize(C %*% x)
  constraints <- list(
    x >= 0,
    x[1] + x[2] <= 5,
    2 * x[1] - x[2] >= 0,
    -x[1] + 3 * x[2] >= 0,
    x[3] + x[4] >= 0.5,
    x[3] >= 1.1
  )
  problem <- Problem(objective, constraints)
  opt_val <- psolve(problem, solver = "HIGHS")
  x_val <- as.numeric(value(x))
  expect_true(is.finite(opt_val))
  ## x3 should be integer
  expect_equal(x_val[3], round(x_val[3]), tolerance = 1e-6)
  ## All nonneg
  expect_true(all(x_val >= -1e-6))
  ## x3 >= 1.1, so integer => x3 >= 2
  expect_true(x_val[3] >= 2 - 1e-6)
})

# ═══════════════════════════════════════════════════════════════════
# 9. sparse-inverse-covariance-estimation.Rmd — SDP (log_det)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("sparse inverse covariance: log_det + PSD variable", {
  skip_if_not_installed("clarabel")
  set.seed(1)
  n <- 5  # smaller than docs (10) for speed
  m <- 200
  A <- matrix(rnorm(n * n) * (runif(n * n) < 0.3), n, n)
  Strue <- A %*% t(A) + 0.05 * diag(n)
  R <- base::solve(Strue)
  ## Sample from multivariate normal
  sqrtR <- chol(R)
  x_sample <- matrix(rnorm(n * m), nrow = m) %*% sqrtR
  Q <- cov(x_sample)

  S <- Variable(c(n, n), PSD = TRUE)
  obj <- Maximize(log_det(S) - matrix_trace(S %*% Q))
  alpha <- 10
  prob <- Problem(obj, list(sum(abs(S)) <= alpha))
  opt_val <- psolve(prob)
  S_val <- value(S)
  expect_true(is.finite(opt_val))
  ## Solution should be symmetric and PSD
  expect_equal(S_val, t(S_val), tolerance = 1e-6)
  eigenvals <- eigen(S_val, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvals > -1e-6))
})

# ═══════════════════════════════════════════════════════════════════
# 10. gentle_intro.Rmd — LP + constraints
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("gentle intro: SVM-style classification", {
  skip_if_not_installed("clarabel")
  ## Simple LP: maximize subject to box constraints
  n <- 2
  x <- Variable(n)
  prob <- Problem(Maximize(sum(x)), list(x >= 0, x <= 3))
  opt_val <- psolve(prob)
  expect_equal(opt_val, 6, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(3, 3), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# 11. solver-parameters.Rmd — Multi-solver
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("same LP solved by all installed solvers", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))

  for (solver_name in installed_solvers()) {
    opt <- psolve(prob, solver = solver_name)
    expect_equal(opt, 2, tolerance = 1e-3,
                 label = paste("solver:", solver_name))
  }
})

# ═══════════════════════════════════════════════════════════════════
# 12. Backward-compat aliases
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("multiply() is deprecated alias for *", {
  rlang::reset_warning_verbosity("cvxr_multiply_deprecated")
  x <- Variable(2)
  expect_warning(e <- multiply(c(2, 3), x), "deprecated")
  expect_true(S7::S7_inherits(e, Expression))
})
