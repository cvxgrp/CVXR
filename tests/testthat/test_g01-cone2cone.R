context("test-g01-cone2cone")
TOL <- 1e-^6

TestDualize.simulate_chain <- function(in_prob) {
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

test_that("test lp 1", {
  # Typical LP.
  sth <- lp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("test lp 2", {
  # Typical LP.
  sth <- lp_2()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("test lp 3", {
  # Unbounded LP.
  sth <- lp_3()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
})

test_that("test lp 4", {
  # Infeasible LP.
  sth <- lp_4()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
})

test_that("test lp 5", {
  # LP with redundant constraints.
  sth <- lp_5()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  check_primal_feasibility(sth, result, tolerance = 1e-4)
  check_complementarity(sth, result, tolerance = 1e-4)
  check_dual_domains(sth, result, tolerance = 1e-4)
})

test_that("test socp 0", {
  sth <- socp_0()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

test_that("test socp 1", {
  sth <- socp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

test_that("test socp 2", {
  sth <- socp_2()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

.socp_3 <- function(axis) {
  sth <- socp_3(axis)
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
}

test_that("test socp3 axis 1", {
  .socp_3(1)
})

test_that("test socp3 axis 2", {
  .socp_3(2)
})

test_that("test expcone 1", {
  sth <- expcone_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("test expcone socp 1", {
  sth <- expcone_socp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("test pcp 2", {
  sth <- pcp_2()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_dual_values(sth, result, tolerance = 1e-3)
})

###############
#             #
# Test Slacks #
#             #
###############

AFF_LP_CASES <- list(c(NONNEG), c())
AFF_SOCP_CASES <- list(c(NONNEG, SOC), c(NONNEG), c(SOC), c())
AFF_EXP_CASES <- list(c(NONNEG, EXP), c(NONNEG), c(EXP), c())
AFF_PCP_CASES <- list(c(NONNEG), c(POW3D), c())
AFF_MIXED_CASES <- list(c(NONNEG), c())

TestSlacks.simulate_chain <- function(in_prob, affine, ...) {
  # Get a ParamConeProg object.
  reductions <- list(Dcp2Cone(), CvxAttr2Constr(), ConeMatrixStuffing())
  chain <- Chain(NULL, reductions)
  tmp <- perform(chain, in_prob)
  cone_prog <- tmp[[1]]
  inv_prob2cone <- tmp[[2]]
  
  # Apply the Slacks reduction, reconstruct a high-level problem, solve the problem, invert the reduction.
  cone_prog <- format_constraints(ConicSolver(), cone_prog, exp_cone_order = c(0, 1, 2))
  # TODO: Finish this.
}

# TODO: Add rest of cone2cone tests and classes like TestSlacks, TestPowND, etc.

