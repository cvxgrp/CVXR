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

test_that("TestDualize: test lp 1", {
  # Typical LP.
  sth <- lp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test lp 2", {
  # Typical LP.
  sth <- lp_2()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test lp 3", {
  # Unbounded LP.
  sth <- lp_3()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test lp 4", {
  # Infeasible LP.
  sth <- lp_4()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test lp 5", {
  # LP with redundant constraints.
  sth <- lp_5()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  check_primal_feasibility(sth, result, tolerance = 1e-4)
  check_complementarity(sth, result, tolerance = 1e-4)
  check_dual_domains(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test socp 0", {
  sth <- socp_0()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test socp 1", {
  sth <- socp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test socp 2", {
  sth <- socp_2()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
})

TestDualize.socp_3 <- function(axis) {
  sth <- socp_3(axis)
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
}

test_that("TestDualize: test socp3 axis 1", {
  TestDualize.socp_3(1)
})

test_that("TestDualize: test socp3 axis 2", {
  TestDualize.socp_3(2)
})

test_that("TestDualize: test expcone 1", {
  sth <- expcone_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test expcone socp 1", {
  sth <- expcone_socp_1()
  result <- TestDualize.simulate_chain(sth@prob)
  verify_objective(sth, result, tolerance = 1e-4)
  verify_primal_values(sth, result, tolerance = 1e-4)
  verify_dual_values(sth, result, tolerance = 1e-4)
})

test_that("TestDualize: test pcp 2", {
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
  tmp <- Slacks.perform(cone_prog, affine)
  data <- tmp[[1]]
  inv_data <- tmp[[2]]
  G <- data[[A_KEY]]
  h <- data[[B_KEY]]
  f <- data[[C_KEY]]
  K_dir <- data$K_dir
  K_aff <- data$K_aff
  
  G <- Matrix(G, sparse = TRUE)
  y <- Variable(ncol(G))
  objective <- Minimize(t(f) %*% y)
  aff_con <- TestSlacks.set_affine_constraints(G, h, y, K_aff)
  dir_con <- TestSlacks.set_direct_constraints(y, K_dir)
  int_con <- TestSlacks.set_integer_constraints(y, data)
  constraints <- c(aff_con, dir_con, int_con)
  slack_prob <- Problem(objective, constraints)
  result <- solve(slack_prob, ...)
  
  slack_prims <- list()
  slack_prims[[Cone2Cone.FREE]] <- result$getValue(y[1:size(cone_prog@x)])
  slack_sol <- Solution(result$status, result$value, slack_prims, NULL, list())
  cone_sol <- Slacks.invert(slack_sol, inv_data)
  
  # Pass solution up the solving chain.
  in_prob_sol <- invert(chain, cone_sol, inv_prob2cone)
  unpack(in_prob, in_prob_sol)
}

TestSlacks.set_affine_constraints(G, h, y, K_aff) {
  constraints <- list()
  i <- 1
  if(Cone2Cone.ZERO %in% names(K_aff)) {
    dim <- K_aff[[Cone2Cone.ZERO]]
    constraints <- c(constraints, G[i:(i+dim-1),] %*% y == h[i:(i+dim-1)])
    i <- i + dim
  }
  
  if(Cone2Cone.NONNEG %in% names(K_aff)) {
    dim <- K_aff[[Cone2Cone.NONNEG]]
    constraints <- c(constraints, G[i:(i+dim-1),] %*% y <= h[i:(i+dim-1)])
    i <- i + dim
  }
  
  for(dim in K_aff[[Cone2Cone.SOC]]) {
    expr <- h[i:(i+dim-1)] - G[i:(i+dim-1),] %*% y
    constraints <- c(constraints, SOC(expr[1], expr[2:nrow(expr)]))
    i <- i + dim
  }
  
  if(Cone2Cone.EXP %in% names(K_aff)) {
    dim <- 3*K_aff[[Cone2Cone.EXP]]
    expr <- reshape_expr(h[i:(i+dim-1)] - G[i:(i+dim-1),] %*% y, c(floor(dim/3), 3), byrow = TRUE)
    constraints <- c(constraints, ExpCone(expr[,1], expr[,2], expr[,3]))
    i <- i + dim
  }
  
  if(Cone2Cone.POW3D %in% names(K_aff)) {
    alpha <- as.vector(K_aff[[Cone2Cone.POW3D]])
    expr <- reshape_expr(h[i:nrow(h)] - G[i:nrow(G),] %*% y, c(length(alpha), 3), byrow = TRUE)
    constraints <- c(constraints, PowCone3D(expr[,1], expr[,2], expr[,3], alpha))
  }
  return(constraints)
}

TestSlacks.set_direct_constraints <- function(y, K_dir) {
  constraints <- list()
  i <- K_dir[[Cone2Cone.FREE]]
  if(Cone2Cone.NONNEG %in% names(K_dir)) {
    dim <- K_dir[[Cone2Cone.NONNEG]]
    constraints <- c(constraints, y[i:(i+dim-1)] >= 0)
    i <- i + dim
  }
  
  for(dim in K_dir[[Cone2Cone.SOC]]) {
    constraints <- c(constraints, SOC(y[i], y[(i+1):(i+dim-1)]))
    i <- i + dim
  }
  
  if(Cone2Cone.EXP %in% names(K_dir)) {
    dim <- 3*K_dir[[Cone2Cone.EXP]]
    expr <- reshape_expr(y[i:(i+dim-1)], c(floor(dim/3), 3), byrow = TRUE)
    constraints <- c(constraints, ExpCone(expr[,1], expr[,2], expr[,3]))
    i <- i + dim
  }
  
  if(Cone2Cone.POW3D %in% names(K_dir)) {
    alpha <- as.vector(K_dir[[Cone2Cone.POW3D]])
    expr <- reshape_expr(y[i:nrow(y)], c(length(alpha), 3), byrow = TRUE)
    constraints <- c(constraints, PowCone3D(expr[,1], expr[,2], expr[,3], alpha))
  }
  return(constraints)
}

TestSlacks.set_integer_constraints <- function(y, data) {
  constraints <- list()
  if(BOOL_IDX %in% names(data)) {
    expr <- y[data[[BOOL_IDX]]]
    z <- Variable(size(expr), boolean = TRUE)
    constraints <- c(constraints, expr == z)
  }
  
  if(data[[INT_IDX]]) {
    expr <- y[data[[INT_IDX]]]
    z <- Variable(size(expr), integer = TRUE)
    constraints <- c(constraints, expr == z)
  }
  return(constraints)
}

test_that("TestSlacks: test lp 2", {
  # Typical LP.
  sth <- lp_2()
  for(affine in AFF_LP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test lp 3", {
  # Unbounded LP.
  sth <- lp_3()
  for(affine in AFF_LP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test lp 4", {
  # Infeasible LP.
  sth <- lp_4()
  for(affine in AFF_LP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test socp 2", {
  sth <- socp_2()
  for(affine in AFF_SOCP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test socp 3", {
  for(axis in c(1, 2)) {
    sth <- socp_3(axis)
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test expcone 1", {
  sth <- expcone_1()
  for(affine in AFF_EXP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test expcone socp 1", {
  sth <- expcone_socp_1()
  for(affine in AFF_MIXED_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

test_that("TestSlacks: test pcp 1", {
  sth <- pcp_1()
  for(affine in AFF_PCP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "SCS", eps = 1e-8)
    verify_objective(sth, result, tolerance = 1e-3)
    verify_primal_values(sth, result, tolerance = 1e-3)
  }
})

test_that("TestSlacks: test pcp 2", {
  sth <- pcp_2()
  for(affine in AFF_PCP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "SCS", eps = 1e-8)
    verify_objective(sth, result, tolerance = 1e-3)
    verify_primal_values(sth, result, tolerance = 1e-3)
  }
})

test_that("TestSlacks: test mi lp 1", {
  sth <- mi_lp_1()
  for(affine in AFF_LP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS_BB")
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_values(sth, result, tolerance = 1e-4)
  }
})

# TODO: Known bug in ECOS BB.
# test_that("TestSlacks: test mi socp 1", {
#  sth <- mi_socp_1()
#  for(affine in AFF_SOCP_CASES) {
#    result <- TestSlacks.simulate_chain(sth@prob, affine, solver = "ECOS_BB")
#    verify_objective(sth, result, tolerance = 1e-4)
#    verify_primal_values(sth, result, tolerance = 1e-4)
#  }
# })

test_that("TestSlacks: test mi socp 2", {
  can_solve <- FALSE
  for(svr in INSTALLED_MI) {
    if(svr %in% MI_SOCP && svr != ECOS_BB) {
      can_solve <- TRUE
      break
    }
  }
  
  if(!can_solve) {
    print("No appropriate mixed-integer SOCP solver is installed")
    return()
  }
  
  sth <- mi_socp_2()
  for(affine in AFF_SOCP_CASES) {
    result <- TestSlacks.simulate_chain(sth@prob, affine)
    verify_objective(sth, result, tolerance = 1e-4)
    verify_primal_value(sth, result, tolerance = 1e-4)
  }
})

##############
#            #
# Test PowND #
#            #
##############

TestPowND.pcp_3 <- function(axis) {
  # A modification of pcp_2. Reformulate
  #
  #         max  (x^0.2)*(y^0.8) + z^0.4 - x
  #          s.t. x + y + z/2 == 2
  #               x, y, z >= 0
  #      Into
  #
  #          max  x3 + x4 - x0
  #          s.t. x0 + x1 + x2 / 2 == 2,
  #
  #               W := [[x0, x2],
  #                     [x1, 1.0]]
  #               z := [x3, x4]
  #               alpha := [[0.2, 0.4],
  #                         [0.8, 0.6]]
  #               (W, z) in PowND(alpha, axis = 2)
  x <- Variable(3)
  expect_x <- c(0.06393515, 0.78320961, 2.30571048)
  hypos <- Variable(2)
  expect_hypos <- NULL
  objective <- Maximize(sum(hypos) - x[1])
  W <- bmat(list(list(x[1], x[3]),
                 list(x[2], 1.0)))
  alpha <- rbind(c(0.2, 0.4), c(0.8, 0.6))
  
  if(axis == 1) {
    W <- t(W)
    alpha <- t(alpha)
  }
  
  con_pairs <- list(list(x[1] + x[2] + 0.5*x[3] == 2, NULL),
                    list(PowConeND(W, hypos, alpha, axis = axis), NULL))
  obj_pairs <- list(objective, 1.8073406786220672)
  var_pairs <- list(list(x, expect_x),
                    list(hypos, expect_hypos))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("TestPowND: test pcp 3a", {
  sth <- TestPowND.pcp_3(axis = 2)
  result <- solve(sth, solver = "SCS", eps = 1e-8)
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_complementarity(sth, result, tolerance = 1e-3)
})

test_that("TestPowND: test pcp 3b", {
  sth <- TestPowND.pcp_3(axis = 1)
  result <- solve(sth, solver = "SCS", eps = 1e-8)
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_complementarity(sth, result, tolerance = 1e-3)
})

TestPowND.pcp_4 <- function(ceei = TRUE) {
  # A power cone formulation of a Fisher market equilibrium pricing model.
  # ceei = Competitive Equilibrium from Equal Incomes.
  
  # Generate test data.
  set.seed(0)
  n_buyer <- 4
  n_items <- 6
  V <- matrix(runif(n_buyer*n_items), nrow = n_buyer, ncol = n_items)
  X <- Variable(n_buyer, n_items, nonneg = TRUE)
  u <- sum(multiply(V, X), axis = 1)
  
  if(ceei)
    b <- rep(1, n_buyer)/n_buyer
  else
    b <- c(0.3, 0.15, 0.2, 0.35)
  
  log_objective <- Maximize(sum(multiply(b, log(u))))
  log_cons <- list(sum(X, axis = 2) <= 1)
  log_prob <- Problem(log_objective, log_cons)
  result <- solve(log_prob, solver = "SCS", eps = 1e-8)
  expect_X <- result$getValue(X)
  
  z <- Variable()
  pow_objective <- list(Maximize(z), exp(result$getValue(log_prob)))
  pow_cons <- list(list(sum(X, axis = 2) <= 1, NULL),
                   list(PowConeND(W = u, z = z, alpha = b), NULL))
  pow_vars <- list(list(X, expect_X))
  sth <- SolverTestHelper(pow_objective, pow_vars, pow_cons)
  return(sth)
}

test_that("TestPowND: test pcp 4a", {
  sth <- TestPowND.pcp_4(ceei = TRUE)
  result <- solve(sth, solver = "SCS", eps = 1e-8)
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_complementarity(sth, result, tolerance = 1e-3)
})

test_that("TestPowND: test pcp 4b", {
  sth <- TestPowND.pcp_4(ceei = FALSE)
  result <- solve(sth, solver = "SCS", eps = 1e-8)
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_complementarity(sth, result, tolerance = 1e-3)
})

######################
#                    #
# Test Rel Entr Quad #
#                    #
######################

TestRelEntrQuad.expcone_1 <- function() {
  # min   3 * x[1] + 2 * x[2] + x[3]
  # s.t.  0.1 <= x[1] + x[2] + x[3] <= 1
  #       x >= 0
  #       and ...
  #         x[1] >= x[2] * exp(x[3] / x[2])
  #       equivalently ...
  #         x[1] / x[2] >= exp(x[3] / x[2])
  #         log(x[1] / x[2]) >= x[3] / x[2]
  #         x[2] log(x[2] / x[1]) <= -x[3]
  
  x <- Variable(3, 1)
  cone_con <- as_quad_approx(ExpCone(x[3], x[2], x[1]), 5, 5)
  constraints <- list(sum(x) <= 1.0, sum(x) >= 0.1, x >= 0, cone_con)
  obj <- Minimize(3*x[1] + 2*x[2] + x[3])
  obj_pair <- list(obj, 0.23534820622420757)
  expect_exp <- list(c(-1.35348213), c(-0.35348211), c(0.64651792))
  con_pairs <- list(list(constraints[[1]], 0),
                    list(constraints[[2]], 2.3534821130067614),
                    list(constraints[[3]], matrix(0, nrow = 3, ncol = 1)),
                    list(constraints[[4]], expect_exp))
  expect_x <- matrix(c(0.05462721, 0.02609378, 0.01927901), nrow = 3, ncol = 1)
  var_pairs <- list(x, expect_x)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("TestRelEntrQuad: expcone 1", {
  sth <- TestRelEntrQuad.expcone_1()
  result <- solve(sth, solver = "ECOS")
  verify_primal_values(sth, result, tolerance = 1e-2)
  verify_objective(sth, result, tolerance = 1e-2)
})

TestRelEntrQuad.expcone_socp_1 <- function() {
  # A random risk-parity portfolio optimization problem.
  sigma <- rbind(c(1.83, 1.79, 3.22),
                 c(1.79, 2.18, 3.18),
                 c(3.22, 3.18, 8.69))
  L <- t(chol(sigma))
  c <- 0.75
  t <- Variable(name = "t")
  x <- Variable(3, name = "x")
  s <- Variable(3, name = "s")
  e <- Constant(matrix(1, nrow = 3, ncol = 1))
  
  objective <- Minimize(t - c * t(e) %*% s)
  con1 <- norm2(t(L) %*% x) <= t
  con2 <- as_quad_approx(ExpCone(s, e, x), 5, 5)
  
  # SolverTestHelper data.
  obj_pair <- list(objective, 4.0751197)
  var_pairs <- list(
    list(x, c(0.57608346, 0.54315695, 0.28037716)),
    list(s, c(-0.55150, -0.61036, -1.27161))
  )
  con_pairs <- list(list(con1, 1.0),
                    list(con2, list(NULL, NULL, NULL)))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("TestRelEntrQuad: expcone socp 1", {
  sth <- TestRelEntrQuad.expcone_socp_1()
  result <- solve(sth, solver = "ECOS")
  verify_primal_values(sth, result, tolerance = 1e-2)
  verify_objective(sth, result, tolerance = 1e-2)
})

sdp_ipm_installed <- function() {
  viable <- intersect(c("CVXOPT", "MOSEK", "COPT"), installed_solvers())
  return(length(viable) > 0)
}

# TODO: Add TestOpRelConeQuad functions from cone2cone unit tests.

