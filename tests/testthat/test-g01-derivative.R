contex("test-g01-derivative")
TOL <- 1e-6

SOLVE_METHODS <- c("SCS", "ECOS")
EPS_NAME = list(SCS = "eps", ECOS = "abstol")

###################################################
#                                                 #
#                Test Backward                    #
#                                                 #
# Test backward(problem) and derivative(problem). #
#                                                 #
###################################################

test_that("test scalar quadratic", {
  b <- Parameter()
  x <- Variable()
  quadratic <- square(x - 2*b)
  problem <- Problem(Minimize(quadratic), list(x >= 0))
  value(b) <- 3
  result <- solve(problem, solver = "DIFFCP", requires_grad = TRUE, eps = 1e-10)
  expect_equal(result$getValue(x), 6, tolerance = TOL)
  back <- backward(problem)
  
  # x* = 2*b, dx*/db = 2.
  # gradient(x) == NA defaults to 1.0.
  expect_equal(back$getGradient(b), 2, tolerance = TOL)
  gradient(x) <- 4
  # quadratic <- square(x - 2*b)
  # problem <- Problem(Minimize(quadratic), list(x >= 0))
  back <- backward(problem)
  expect_equal(back$getGradient(b), 8, tolerance = TOL)
  gradcheck(problem, atol = 1e-4)
  perturbcheck(problem, atol = 1e-4)
  
  result <- solve(problem, solver = "DIFFCP", requires_grad = TRUE, eps = 1e-10)
  delta(b) <- 1e-3
  deriv <- derivative(problem)
  expect_equal(deriv$getDelta(x), 2e-3, tolerance = TOL)
})

test_that("test l1 square", {
  set.seed(0)
  n <- 3
  x <- Variable(n)
  A <- Parameter(n, n)
  b <- Parameter(n, name = "b")
  objective <- Minimize(p_norm(A %*% x - b, p = 1))
  problem <- Problem(objective)
  expect_true(is_dpp(problem))
  
  L <- matrix(rnorm(n*n), nrow = n, ncol = n)
  value(A) <- t(L) %*% L + diag(n)
  value(b) <- matrix(rnorm(n))
  gradcheck(problem)
  perturbcheck(problem)
})

test_that("test l1 rectangle", {
  set.seed(0)
  m <- 3
  n <- 2
  x <- Variable(n)
  A <- Parameter(m, n)
  b <- Parameter(m, name = "b")
  objective <- Minimize(p_norm(A %*% x - b, p = 1))
  problem <- Problem(objective)
  expect_true(is_dpp(problem))
  
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  gradcheck(problem, atol = 1e-3)
  perturbcheck(problem, atol = 1e-3)
})

test_that("test least squares", {
  set.seed(0)
  m <- 20
  n <- 5
  A <- Parameter(m, n)
  b <- Parameter(m)
  x <- Variable(n)
  obj <- sum_squares(A %*% x - b) + sum_squares(x)
  problem <- Problem(Minimize(obj))
  
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  gradcheck(problem, solve_methods = c("SCS"))
  perturbcheck(problem, solve_methods = c("SCS"))
})

test_that("test logistic regression", {
  set.seed(0)
  N <- 5
  n <- 2
  X_np <- matrix(rnorm(N*n), nrow = N, ncol = n)
  a_true <- matrix(rnorm(n), nrow = n, ncol = 1)
  
  sigmoid <- function(z) {
    return(1 / (1 + exp(-z)))
  }
  
  y <- round(sigmoid(X_np %*% a_true + rnorm(N)*0.5))
  
  a <- Variable(n, 1)
  X <- Parameter(N, n)
  lam <- Parameter(nonneg = TRUE)
  log_likelihood <- sum(multiply(y, X %*% a) - t(log_sum_exp(t(hstack(matrix(0, nrow = N, ncol = 1), X %*% a)), axis = 2, keepdims = TRUE)))
  problem <- Problem(Minimize(-log_likelihood + lam*sum_squares(a)))
  value(X) <- X_np
  value(lam) <- 1
  
  # TODO: Too low, but this problem is ill-conditioned.
  gradcheck(problem, solve_methods = c("SCS"), atol = 1e-1, eps = 1e-8)
  perturbcheck(problem, solve_methods = c("SCS"), atol = 1e-4)
})

test_that("test entropy maximization", {
  set.seed(0)
  n <- 5
  m <- 3
  p <- 2
  
  tmp <- matrix(runif(n))
  A_np <- matrix(rnorm(m*n), nrow =  m, ncol = n)
  b_np <- A_np %*% tmp
  F_np <- matrix(rnorm(p, n), nrow = p, ncol = n)
  g_np <- F_np %*% tmp + runif(p)
  
  x <- Variable(n)
  A <- Parameter(m, n)
  b <- Parameter(m)
  Fp <- Parameter(p, n)
  g <- Parameter(p)
  obj <- Maximize(sum(entr(x)) - sum_squares(x))
  constraints <- list(A %*% x == b, Fp %*% x <= g)
  problem <- Problem(obj, constraints)
  
  value(A) <- A_np
  value(b) <- b_np
  value(Fp) <- F_np
  value(g) <- g_np
  
  gradcheck(problem, solve_methods = c("SCS"), atol = 1e-2, eps = 1e-8, max_iters = 1e4)
  perturbcheck(problem, solve_methods = c("SCS"), atol = 1e-4)
})

test_that("test lml", {
  set.seed(0)
  k <- 2
  x <- Parameter(4)
  y <- Variable(4)
  obj <- -t(x) %*% y - sum(entr(y)) - sum(entr(1 - y))
  cons <- list(sum(y) == k)
  problem <- Problem(Minimize(obj), cons)
  
  value(x) <- c(1, -1, -1, -1)
  
  # TODO: This tolerance is too low.
  gradcheck(problem, solve_methods = c("SCS"), atol = 1e-2)
  perturbcheck(problem, solve_methods = c("SCS"), atol = 1e-4)
})

test_that("test sdp", {
  set.seed(0)
  n <- 3
  p <- 3
  C <- Parameter(n, n)
  As <- lapply(1:p, function(i) { Parameter(n, n) })
  bs <- lapply(1:p, function(i) { Parameter(1, 1) })
  
  value(C) <- matrix(rnorm(n*n), nrow = n, ncol = n)
  for(i in 1:p) {
    value(As[[i]]) <- matrix(rnorm(n*n), nrow = n, ncol = n)
    value(bs[[i]]) <- matrix(rnorm(1))
  }
  
  X <- Variable(n, n, PSD = TRUE)
  constraints <- lapply(1:p, function(i) { matrix_trace(As[[i]] %*% X) == bs[[i]] })
  problem <- Problem(Minimize(matrix_trace(C %*% X) + sum_squares(X)), constraints)
  gradcheck(problem, solve_methods = c("SCS"), atol = 1e-3, eps = 1e-10)
  perturbcheck(problem, solve_methods = c("SCS"))
})

test_that("test forget requires grad", {
  set.seed(0)
  m <- 20
  n <- 5
  A <- Parameter(m, n)
  b <- Parameter(m)
  x <- Variable(n)
  
  obj <- sum_squares(A %*% x - b) + sum_squares(x)
  problem <- Problem(Minimize(obj))
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  result <- solve(problem, solver = "SCS")
  
  expect_error(backward(problem), "backward can only be called after calling solve with requires_grad = TRUE")
  expect_error(derivative(problem), "derivative can only be called after calling solve with requires_grad = TRUE")
})

test_that("test infeasible", {
  x <- Variable()
  param <- Parameter()
  problem <- Problem(Minimize(param), list(x >= 1, x <= -1))
  value(param) <- 1
  result <- solve(problem, solver = "DIFFCP", requires_grad = TRUE)
  
  expect_error(backward(problem), "Backpropagating through infeasible/unbounded")
  expect_error(derivative(problem), "Differentiating through infeasible/unbounded")
})

test_that("test unbounded", {
  x <- Variable()
  param <- Parameter()
  problem <- Problem(Minimize(x), list(x <= param))
  value(param) <- 1
  result <- solve(problem, solver = "DIFFCP", requires_grad = TRUE)
  
  expect_error(backward(problem), "Backpropagating through infeasible/unbounded")
  expect_error(derivative(problem), "Differentiating through infeasible/unbounded")
})

test_that("test unsupported solver", {
  x <- Variable()
  param <- Parameter()
  problem <- Problem(Minimize(x), list(x <= param))
  value(param) <- 1
  expect_error(solve(problem, solver = "ECOS", requires_grad = TRUE), 
               "When requires_grad = TRUE, the only supported solver is SCS")
})

test_that("test zero in problem data", {
  x <- Variable()
  param <- Parameter()
  value(param) <- 0.0
  problem <- Problem(Minimize(x), list(param*x >= 0))
  data <- get_problem_data(problem, "DIFFCP")[[1]]
  A <- data[[A_KEY]]
  expect_in(0.0, A$data)
})

