contex("test-g01-derivative")
TOL <- 1e-6

SOLVE_METHODS <- c("SCS", "ECOS")
EPS_NAME = list(SCS = "eps", ECOS = "abstol")

perturbcheck <- function(problem, gp = FALSE, solve_methods = SOLVE_METHODS, delta = 1e-5, atol = 1e-6, eps = 1e-9, ...) {
  # Checks the analytical derivative against a numerical computation.
  for(solver in solve_methods) {
    set.seed(0)
    eps_opt <- list()
    eps_opt[[EPS_NAME[[solver]]]] <- eps
    if(length(parameters(problem)) == 0) {
      result <- do.call("solve", c(list(a = problem, solver = "DIFFCP", gp = gp, requires_grad = TRUE, solve_method = solver),
                                   eps_opt, list(...)))
      deriv <- derivative(problem)
      for(variable in variables(problem))
        expect_equal(deriv$getDelta(variable), 0.0)
    }
    
    # Compute perturbations analytically.
    for(param in parameters(problem))
      delta(param) <- delta*matrix(rnorm(size(param)), nrow = nrow(param), ncol = ncol(param))
    result <- do.call("solve", c(list(a = problem, solver = "DIFFCP", gp = gp, requires_grad = TRUE, solve_method = solver),
                                 eps_opt, list(...)))
    deriv <- derivative(problem)
    variable_values <- lapply(variables(problem), value)
    deltas <- lapply(variables(problem), delta)
    
    # Compute perturbations numerically.
    old_values <- list()
    for(param in parameters(problem)) {
      old_values[[as.character(id(param))]] <- value(param)
      value(param) <- value(param) + delta(param)
    }
    result <- do.call("solve", c(list(a = problem, solver = solver, gp = gp), eps_opt, list(...)))
    prob_vars <- variables(problem)
    num_deltas <- lapply(1:length(variable_values), function(i) {
      v <- prob_vars[[i]]
      old_value <- variable_values[[i]]
      value(v) - old_value
    })
    
    for(i in 1:length(deltas)) {
      analytical <- deltas[[i]]
      numerical <- num_deltas[[i]]
      expect_true(is.allclose(analytical, numerical, atol = atol))
    }
    
    for(param in parameters(problem))
      value(param) <- old_values[[as.character(id(param))]]
  }
}

gradcheck <- function(problem, gp = FALSE, solve_methods = SOLVE_METHODS, delta = 1e-5, atol = 1e-4, eps = 1e-9, ...) {
  # Checks the analytical adjoint derivative against a numerical computation.
  for(solver in solve_methods) {
    eps_opt <- list()
    eps_opt[[EPS_NAME[[solver]]]] <- eps
    # Default of 15k iterations for SCS.
    if(solver == "SCS" && !("max_iters" %in% names(list(...))))
      kwargs$max_iters <- 15000
    
    size <- sum(sapply(parameters(problem), size))
    values <- rep(0, size)
    offset <- 0
    
    for(param in parameters(problem)) {
      values[(offset + 1):(offset + size(param) + 1)] <- as.vector(value(param))
      value(param) <- matrix(values[(offset + 1):(offset + size(param) + 1)], nrow = nrow(param), ncol = ncol(param), byrow = TRUE)
      offset <- offset + size(param)
    }
    
    numgrad <- matrix(0, nrow = nrow(values), ncol = ncol(values))
    for(i in 1:size(values)) {
      old <- values[i]
      values[i] <- old + 0.5*delta
      result <- do.call("solve", c(list(a = problem, solver = solver, gp = gp), eps_opt, list(...)))
      left_solns <- do.call("c", lapply(variables(problem), function(x) { as.vector(result$getValue(x)) }))
      
      values[i] <- old - 0.5*delta
      result <- do.call("solve", c(list(a = problem, solver = solver, gp = gp), eps_opt, list(...)))
      right_solns <- do.call("c", lapply(variables(problem), function(x) { as.vector(result$getValue(x)) }))
      
      num_grad[i] <- (sum(left_solns) - sum(right_solns))/delta
      values[i] <- old
    }
    
    numgrads <- list()
    offset <- 0
    for(param in parameters(problem)) {
      numgrads <- c(numgrads, matrix(numgrad[(offset + 1):(offset + size(param) + 1)], nrow = nrow(param), ncol = ncol(param), byrow = TRUE))
      offset <- offset + size(param)
    }
    
    old_gradients <- list()
    for(x in variables(problem)) {
      old_gradients[[as.character(id(x))]] <- gradient(x)
      gradient(x) <- NA_real_
    }
    result <- do.call("solve", c(list(a = problem, solver = "DIFFCP", requires_grad = TRUE, gp = gp, solve_method = solver), eps_opt, list(...)))
    back <- backward(problem)
    
    prob_parms <- parameters(problem)
    for(i in 1:length(numgrads)) {
      param <- prob_parms[[i]]
      numgrad <- numgrads[[i]]
      expect_true(is.allclose(back$getGradient(param), numgrad, atol = atol))
    }
    
    for(x in variables(problem))
      gradient(x) <- old_gradients[[as.character(id(x))]]
  }
}

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

###################################################
#                                                 #
#              Test Backward DGP                  #
#                                                 #
# Test backward(problem) and derivative(problem). #
#                                                 #
###################################################

test_that("test one minus analytic", {
  # Construct a problem with solution
  # x^\star(\alpha) = 1 - \alpha^2, and derivative
  # x^\star'(\alpha) = -2\alpha
  
  alpha <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE)
  objective <- Maximize(x)
  constr <- list(one_minus_pos(x) >= alpha^2)
  problem <- Problem(objective, constr)
  
  value(alpha) <- 0.4
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-5)
  expect_equal(result$getValue(x), 1 - 0.4^2, tolerance = 1e-3)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(back$getGradient(alpha), -2*0.4, tolerance = 1e-3)
  expect_equal(deriv$getDelta(x), -2*0.4*1e-5, tolerance = 1e-3)
  
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  
  value(alpha) <- 0.5
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-5)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(result$getValue(x), 1 - 0.5^2, tolerance = 1e-3)
  expect_equal(back$getGradient(alpha), -2*0.5, tolerance = 1e-3)
  expect_equal(deriv$getDelta(x), -2*0.5*1e-5, tolerance = 1e-3)
  
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = tRUE, solve_methods = c("SCS"), atol = 1e-3)
})

test_that("test analytic param in exponent", {
  # Construct a problem with solution
  # x^\star(\alpha) = 1 - \alpha^2, and derivative
  # x^\star'(\alpha) = -log(2)*2^\alpha
  
  base <- 2.0
  alpha <- Parameter()
  x <- Variable(pos = TRUE)
  objective <- Maximize(x)
  constr <- list(one_minus_pos(x) >= Constant(base)^alpha)
  problem <- Problem(objective, constr)
  
  value(alpha) <- -1.0
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-6)
  expect_equal(result$getValue(x), 1 - base^(-1.0), tolerance = TOL)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(back$getGradient(alpha), -log(base)*base^(-1.0), tolerance = TOL)
  expect_equal(deriv$getDelta(x), back$getGradient(alpha)*1e-5, tolerance = 1e-3)
  
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  
  value(alpha) <- -1.2
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-6)
  expect_equal(result$getValue(x), 1 - base^(-1.2), tolerance = TOL)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(back$getGradient(alpha), -log(base)*base^(-1.2), tolerance = TOL)
  expect_equal(deriv$getDelta(x), back$getGradient(alpha)*1e-5, tolerance = 1e-3)
  
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
})

test_that("test param used twice", {
  # Construct a problem with solution
  # x^\star(\alpha) = 1 - \alpha^2 - \alpha^3, and derivative
  # x^\star'(\alpha) = -2\alpha - 3\alpha^2
  
  alpha <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE)
  objective <- Maximize(x)
  constr <- list(one_minus_pos(x) >= alpha^2 + alpha^3)
  problem <- Problem(objective, constr)
  
  value(alpha) <- 0.4
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-6)
  expect_equal(result$getValue(x), 1 - 0.4^2 - 0.4^3, tolerance = TOL)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(back$getGradient(alpha), -2*0.4 - 3*0.4^2, tolerance = TOL)
  expect_equal(deriv$getDelta(x), back$getGradient(alpha)*1e-5, tolerance = TOL)
  
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
})

test_that("test param used in exponent and elsewhere", {
  # Construct a problem with solution
  # x^\star(\alpha) = 1 - 0.3^\alpha - \alpha^2, and derivative
  # x^\star'(\alpha) = -log(0.3) * 0.2^\alpha - 2*alpha
  
  base <- 0.3
  alpha <- Parameter(pos = TRUE, value = 0.5)
  x <- Variable(pos = TRUE)
  objective <- Maximize(x)
  constr <- list(one_minus_pos(x) >= Constant(base)^alpha + alpha^2)
  problem <- Problem(objective, constr)
  
  delta(alpha) <- 1e-5
  result <- solve(problem, solver = "DIFFCP", gp = TRUE, requires_grad = TRUE, eps = 1e-5)
  expect_equal(result$getValue(x), 1 - base^(0.5) - 0.5^2)
  back <- backward(problem)
  deriv <- derivative(problem)
  expect_equal(back$getGradient(alpha), -log(base)*base^(0.5) - 2*0.5, tolerance = TOL)
  expect_equal(deriv$getDelta(x), back$getGradient(alpha)*1e-5, tolerance = 1e-3, tolerance = TOL)
})

test_that("test basic gp", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  a <- Parameter(pos = TRUE)
  b <- Parameter(pos = TRUE)
  c <- Parameter()
  
  constraints <- list(a*(x*y + x*z + y*z) <= b, x >= y^c)
  problem <- Problem(Minimize(1/(x*y*z)), constraints)
  expect_true(is_dgp(problem, dpp = TRUE))
  
  value(a) <- 2.0
  value(b) <- 1.0
  value(c) <- 0.5
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
})

test_that("test max_elemwise", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  a <- Parameter(value = 0.5)
  b <- Parameter(pos = TRUE, value = 3.0)
  c <- Parameter(pos = TRUE, value = 1.0)
  d <- Parameter(pos = TRUe, value = 4.0)
  
  prod1 <- x*y^a
  prod2 <- b*x*y^a
  obj <- Minimize(max_elemwise(prod1, prod2))
  constr <- list(x == c, y == d)
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE, atol = 1e-3)
  perturbcheck(problem, gp = TRUE, atol = 1e-3)
})

test_that("test div", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 3)
  b <- Parameter(pos = TRUE, value = 1)
  problem <- Problem(Minimize(x*y), list(y/a <= x, y >= b))
  gradcheck(problem, gp = TRUE)
  perturbcheck(problem, gp = TRUE)
})

test_that("test one_minus_pos", {
  x <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 3)
  b <- Parameter(pos = TRUE, value = 0.1)
  obj <- Maximize(x)
  constr <- list(one_minus_pos(a*x) >= a*b)
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
})

test_that("test paper example one_minus_pos", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 2)
  b <- Parameter(pos = TRUE, value = 1)
  c <- Parameter(pos = TRUE, value = 3)
  obj <- Minimize(x*y)
  const <- list((y * one_minus_pos(x/y))^a >= b, x >= y/c)
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-3)
  perturbcheck(problem, solve_methods = c("SCS"), gp = TRUE, atol = 1e-3)
})

test_that("test matrix constraint", {
  X <- Variable(2, 2, pos = TRUE)
  a <- Parameter(pos = TRUE, value = 0.1)
  obj <- Minimize(geo_mean(vec(X)))
  constr <- list(diag(X) == a, hstack(X[1,2], X[2,1]) == 2*a)
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE)
  perturbcheck(problem, gp = TRUE)
})

test_that("test paper example exp log", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 0.2)
  b <- Parameter(pos = TRUe, value = 0.3)
  obj <- Minimize(x*y)
  constr <- list(exp(a*y/x) <= log(b*y))
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2, max_iters = 1e4)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2, max_iters = 5000)
})

test_that("test matrix completion", {
  X <- Variable(3, 3, pos = TRUE)
  # TODO: pf matrix completion not differentiable...? I could believe that... or a bug?
  obj <- Minimize(sum(X))
  known_indices <- cbind(c(1, 1, 2, 3, 3), c(1, 3, 2, 1, 2))
  known_values <- c(1.0, 1.9, 0.8, 3.2, 5.9)
  param <- Parameter(length(known_values), pos = TRUE, value = known_values)
  beta <- Parameter(pos = TRUE, value = 1.0)
  constr <- list(X[known_indices] == param,
                 X[1,2]*X[2,1]*X[2,3]*X[3,3] == beta)
  problem <- Problem(obj, constr)
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-4)
})

test_that("test rank one nmf", {
  X <- Variable(3, 3, pos = TRUE)
  x <- Variable(3, pos = TRUE)
  y <- Variable(3, pos = TRUE)
  xy <- vstack(x[1]*y, x[2]*y, x[3]*y)
  a <- Parameter(value = -1.0)
  b <- Parameter(6, pos = TRUE, value = c(1.0, 1.9, 0.8, 3.2, 5.9, 1.0))
  
  R <- max_elemwise(multiply(X, xy^a), multiply(X^a, xy))
  objective <- sum(R)
  constraints <- list(X[1,1] == b[1],
                      X[1,3] == b[2],
                      X[2,2] == b[3],
                      X[3,1] == b[4],
                      X[3,2] == b[5],
                      x[1]*x[2]*x[3] == b[6])
  problem <- Problem(Minimize(objective), constraints)
  # SCS struggles to solve this problem (solved/inaccurate, unless max_iters is very high like 10000).
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2, max_iters = 1000)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2, max_iters = 1000)
})

test_that("test documentation prob", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 4.0)
  b <- Parameter(pos = TRUE, value = 2.0)
  c <- Parameter(pos = TRUE, value = 10.0)
  d <- Parameter(pos = TRUE, value = 1.0)
  
  objective_fn <- x*y*z
  constraints <- list(a*x*y*z + b*x*z <= c, x <= b*y, y <= b*x, z >= d)
  problem <- Problem(Maximize(objective_fn), constraints)
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-2)
})

test_that("test sum_squares vector", {
  alpha <- Parameter(2, pos = TRUE, value = c(1.0, 1.0))
  beta <- Parameter(pos = TRUE, value = 20)
  kappa <- Parameter(pos = TRUE, value = 10)
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  problem <- Problem(Minimize(sum_squares(alpha + h)), 
                     list(multiply(w, h) >= beta, sum(alpha + w) <= kappa))
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-1, max_iters = 1000)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-1, max_iters = 1000)
})

test_that("test sum matrix", {
  w <- Variable(2, 2, pos = TRUE)
  h <- Variable(2, 2, pos = TRUE)
  alpha <- Parameter(pos = TRUE, value = 1.0)
  beta <- Parameter(pos = TRUE, value = 20)
  kappa <- Parameter(pos = TRUE, value = 10)
  problem <- Problem(Minimize(alpha*sum(h)), 
                     list(multiply(w, h) >= beta, sum(w) <= kappa))
  gradcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-1)
  perturbcheck(problem, gp = TRUE, solve_methods = c("SCS"), atol = 1e-1)
})
