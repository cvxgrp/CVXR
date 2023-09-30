context("test-g01-cone2cone")
TOL <- 1e-^6

simulate_chain <- function(in_prob) {
  # Get a ParamConeProg object.
  reduction <- list(Dcp2Cone(), CvxAttr2Constr(), ConeMatrixStuffing())
  chain <- Chain(NULL, reductions)
  tmp <- perform(chain, in_prob)
  cone_prog <- tmp[[1]]
  inv_prob2cone <- tmp[[2]]
  
  # Dualize the problem, reconstruct a high-level CVXR problem for the dual.
  # Solve the problem, invert the dualize reduction.
  cone_prog <- format_constraints(ConicSolver(), cone_prog, exp_cone_order = c(1,2,3))
  tmp <- Dualize.perform(cone_prog)
  data <- tmp[[1]]
  inv_data <- tmp[[2]]
  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  c <- data[[C_KEY]]
  K_dir <- data$K_dir
  
  y <- Variable(ncol(A))
  constraints <- list(A %*% y == b)
  i <- K_dir[[Cone2Cone.FREE]]
  dual_prims <- list()
  dual_prims[[Cone2Cone.FREE]] <- y[seq(i)]
  dual_prims[[Cone2Cone.SOC]] <- list()
  
  if(K_dir[[Cone2Cone.NONNEG]]) {
    dim <- K_dir[[Cone2Cone.NONNEG]]
    dual_prims[[Cone2Cone.NONNEG]] <- y[seq(i+1, i+dim)]
    constraints <- c(constraints, y[seq(i+1, i+dim)] >= 0)
    i <- i + dim
  }
  
  for(dim in K_dir[[Cone2Cone.SOC]]) {
    dual_prims[[Cone2Cone.SOC]] <- c(dual_prims[[Cone2Cone.SOC]], y[seq(i+1, i+dim)])
    if(dim <= 1)
      constraints <- c(constraints, SOC(y[i+1], c()))
    else
      constraints <- c(constraints, SOC(y[i+1], y[seq(i+2, i+dim)]))
    i <- i + dim
  }
  
  if(K_dir[[Cone2Cone.DUAL_EXP]]) {
    exp_len <- 3*K_dir[[Cone2Cone.DUAL_EXP]]
    dual_prims[[Cone2Cone.DUAL_EXP]] <- y[seq(i+1, i+exp_len)]
    y_de <- reshape_expr(y[seq(i+1, i+exp_len)], new_dim = c(floor(exp_len/3), 3), byrow = TRUE)   # Fill rows first.
    constraints <- c(constraints, ExpCone(-y_de[,2], -y_de[,1], exp(1)*y_de[,3]))
    i <- i + exp_len
  }
  
  if(K_dir[[Cone2Cone.DUAL_POW3D]]) {
    alpha <- as.vector(K_dir[[Cone2Cone.DUAL_POW3D]])
    dual_prims[[Cone2Cone.DUAL_POW3D]] <- y[seq(i+1, nrow(y))]
    y_dp <- reshape_expr(y[seq(i+1, nrow(y))], new_dim = c(size(alpha), 3), byrow = TRUE)   # Fill rows first.
    pow_con <- PowCone3D(y_dp[,1]/alpha, y_dp[,2]/(1-alpha), y_dp[,3], alpha)
    constraints <- c(constraints, pow_con)
  }
  
  objective <- Maximize(t(c) %*% y)
  dual_prob <- Problem(objective, constraints)
  result <- solve(dual_prob, solver = "SCS", eps = 1e-8)
  dual_prims[[Cone2Cone.FREE]] <- result$getValue(dual_prims[[Cone2Cone.FREE]])
  
  if(K_dir[[Cone2Cone.NONNEG]])
    dual_prims[[Cone2Cone.NONNEG]] <- result$getValue(dual_prims[[Cone2Cone.NONNEG]])
  dual_prims[[Cone2Cone.SOC]] <- lapply(dual_prims[[Cone2Cone.SOC]], result$getValue)
  
  if(K_dir[[Cone2Cone.DUAL_EXP]])
    dual_prims[[Cone2Cone.DUAL_EXP]] <- result$getValue(dual_prims[[Cone2Cone.DUAL_EXP]])
  
  if(K_dir[[Cone2Cone.DUAL_POW3D]])
    dual_prims[[Cone2Cone.DUAL_POW3D]] <- result$getValue(dual_prims[[Cone2Cone.DUAL_POW3D]])
  
  dual_duals <- list()
  dual_duals[[EQ_DUAL]] <- result$getDualValue(constraints[[1]])
  dual_sol <- Solution(result$status, result$value, dual_prims, dual_duals, list())
  cone_sol <- Dualize.invert(dual_sol, inv_data)
  
  # Pass the solution back up the solving chain.
  in_prob_sol <- invert(chain, cone_sol, inv_prob2cone)
  unpack(in_prob, in_prob_sol)
}
