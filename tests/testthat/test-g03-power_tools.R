context("test-g03-power_tools")

test_that("test multi step dyad completion", {
  # Consider four market equilibrium problems.
  # 
  # The budgets "b" in these problems are chosen so that canonicalization
  # of geo_mean(u, b) hits a recursive code-path in power_tools.dyad_completion(...).
  # 
  # The reference solution is computed by taking the log of the geo_mean objective,
  # which has the effect of making the problem ExpCone representable.
  
  if("MOSEK" %in% installed_solvers())
    log_solve_args <- list(solver = "MOSEK")
  else
    log_solve_args <- list(solver = "ECOS")
  
  n_buyer <- 5
  n_items <- 7
  
  set.seed(0)
  V <- 0.5*(1 + matrix(runif(n_buyer*n_items, nrow = n_buyer, ncol = n_items)))
  X <- Variable(n_buyer, n_items, nonneg = TRUE)
  cons <- list(sum(X, axis = 2) <= 1)
  u <- sum(multiply(V, X), axis = 1)
  bs <- rbind(c(110, 14, 6, 77, 108),
              c(15, 4, 8, 0, 9),
              c(14, 21, 217, 57, 6),
              c(3, 36, 77, 8, 8))
  
  for(i in nrow(bs)) {
    b <- matrix(bs[i,], ncol = 1)
    log_objective <- Maximize(t(b) %*% log(u))
    log_prob <- Problem(log_objective, cons)
    result <- do.call("solve", c(list(a = log_prob), log_solve_args))
    expect_X <- result$getValue(X)
    
    geo_objective <- Maximize(geo_mean(u, b))
    geo_prob <- Problem(geo_objective, cons)
    result_geo <- solve(geo_prob)
    actual_X <- result_geo$getValue(X)
    
    tryCatch({
      expect_equal(actual_X, expect_X, tolerance = 1e-3)
    }, error = function(e) {
      print(paste("Failure at index ", i, " (when b = ", as.character(b), ").", sep = ""))
      result_verb <- do.call("solve", c(list(a = log_prob, verbose = TRUE), log_solve_args))
      print(result_verb$getValue(X))
      result_verb_geo <- solve(geo_prob, verbose = TRUE)
      print(result_verb_geo$getValue(X))
      print("The valuation matrix was")
      print(V)
      stop(e)
    })
  }
})

test_that("test 3d power cone approx", {
  # Use
  #    geo_mean((x,y), (alpha, 1-alpha)) >= |z|
  # as a reformulation of
  #    PowCone3D(x, y, z, alpha).
  # 
  # Check validity of the reformulation by solving
  # orthogonal projection problems.
  
  if("MOSEK" %in% installed_solvers())
    proj_solve_args <- list(solver = "MOSEK")
  else
    proj_solve_args <- list(solver = "SCS", eps = 1e-10)
  
  min_numerator <- 2
  denominator <- 25
  x <- Variable(3)
  
  set.seed(0)
  y <- 10*matrix(runif(3), nrow = 3, ncol = 1)   # The third value doesn't matter.
  
  num_vec <- seq(min_numerator, denominator, 3)
  for(i in seq_along(num_vec)) {
    numerator <- num_vec[i]
    alpha_float <- numerator / denominator
    y[3] <- (y[1]^alpha_float) * (y[2]^(1 - alpha_float)) + 0.05
    objective <- Minimize(p_norm(y - x, 2))
    
    actual_constraints <- list(PowCone3D(x[1], x[2], x[3], alpha_float))
    actual_prob <- Problem(objective, actual_constraints)
    result <- do.call("solve", c(list(a = actual_prob), proj_solve_args))
    actual_x <- result$getValue(x)
    
    weights <- c(alpha_float, 1 - alpha_float)
    approx_constraints <- list(geo_mean(x[1:2], weights) >= abs(x[3]))
    approx_prob <- Problem(objective, approx_constraints)
    result_apx <- solve(approx_prob)
    approx_x <- result_apx$getValue(x)
    
    tryCatch({
      expect_equal(actual_x, approx_x, tolerance = 1e-4)
    }, error = function(e) {
      print(paste("Failure at index ", i, " (when alpha = ", alpha_float, ")", sep = ""))
      stop(e)
    })
  }
})
