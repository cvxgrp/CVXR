lp_0 <- function() {
  x <- Variable(nrow = 2)
  con_pairs <- list(list(x == 0, NULL))
  obj_pairs <- list(Minimize(norm(x, 1) + 1.0), 1)
  var_pairs <- list(list(x, matrix(c(0, 0))))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

lp_1 <- function() {
  x <- Variable(nrow = 2, name = "x")
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= 3,
                      x[1] + 2*x[2] <= 3,
                      x[1] >= 0,
                      x[2] >= 0)
  con_pairs <- list(list(constraints[[1]], 1),
                    list(constraints[[2]], 2),
                    list(constraints[[3]], 3),
                    list(constraints[[3]], 4))
  var_pairs <- list(list(x, matrix(c(1,1))))
  obj_pair <- list(objective, -9)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

lp_2 <- function() {
  x <- Variable(nrow = 2, name = "x")
  objective <- Minimize(x[1] + 0.5*x[2])
  constraints <- list(x[1] >= -100, x[1] <= -10, x[2] == 1)
  con_pairs <- list(list(constraints[[1]], 1),
                    list(constraints[[2]], 0),
                    list(constraints[[3]], -0.5))
  var_pairs <- list(list(x, matrix(c(-100, 1))))
  obj_pair <- list(objective, -99.5)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

lp_3 <- function() {
  # An unbounded problem.
  x <- Variable(5)
  objective <- list(Minimize(sum(x)), -Inf)
  var_pairs <- list(list(x, NULL))
  con_pairs <- list(list(x <= 1, NULL))
  sth <- SolverTestHelper(objective, var_pairs, con_pairs)
  return(sth)
}

lp_4 <- function() {
  # An infeasible problem.
  x <- Variable(5)
  objective <- list(Minimize(sum(x)), Inf)
  var_pairs <- list(list(x, NULL))
  con_pairs <- list(list(x <= 0, NULL), 
                    list(x >= 1, NULL))
  sth <- SolverTestHelper(objective, var_pairs, con_pairs)
  return(sth)
}

lp_5 <- function() {
  # A problem with redundant equality constraints.
  # 10 variables, 6 equality constraints A %*% x == b (two redundant).
  x0 <- matrix(c(0, 1, 0, 2, 0, 4, 0, 5, 6, 7))
  mu0 <- matrix(c(-2, -1, 0, 1, 2, 3.5))
  
  set.seed(0)
  A_min <- matrix(rnorm(4*10), nrow = 4, ncol = 10)
  A_red <- t(A_min) %*% matrix(runif(4*2), nrow = 4, ncol = 2)
  A_red <- t(A_red)
  A <- rbind(A_min, A_red)
  
  b <- A %*% x0   # x0 is primal feasible.
  c_vec <- t(A) %*% mu0   # mu0 is dual feasible.
  c_vec[c(1, 3, 5, 7)] <- c_vec[c(1, 3, 5, 7)] + runif(4)
  
  # c >= t(A) %*% mu0 exhibits complementary slackness with respect to x0.
  # Therefore, (x0, mu0) are primal-dual optimal for...
  x <- Variable(10)
  objective <- list(Minimize(t(c_vec) %*% x), c_vec %*% x0)
  var_pairs <- list(list(x, x0))
  con_pairs <- list(list(x >= 0, NULL), list(A %*% x == b, NULL))
  sth <- SolverTestHelper(objective, var_pairs, con_pairs)
  return(sth)
}

lp_6 <- function() {
  # Test LP with no constraints.
  x <- Variable()
  objective <- Maximize(Constant(0.23)*x)
  obj_pair <- list(objective, Inf)
  var_pairs <- list(list(x, NULL))
  sth <- SolverTestHelper(obj_pair, var_pairs, list())
  return(sth)
}

qp_0 <- function() {
  # Univariate feasible problem.
  x <- Variable(1)
  objective <- Minimize(square(x))
  constraints <- list(x[1] >= 1)
  con_pairs <- list(list(constraints[[1]], 2))
  obj_pair <- list(objective, 1)
  var_pairs <- list(list(x, 1))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

socp_0 <- function() {
  x <- Variable(nrow = 2)
  obj_pair <- list(Minimize(norm(x, 2) + 1), 1)
  con_pairs <- list(list(x == 0, NULL))
  var_pairs <- list(list(x, matrix(c(0,0))))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

socp_1 <- function() {
  # min 3*x[1] + 2*x[2] + x[3]
  # s.t. norm(x,2) <= y,
  #      x[1] + x[2] + 3*x[3] >= 1.0,
  #      y <= 5.
  x <- Variable(nrow = 3)
  y <- Variable()
  soc <- SOC(y, x)
  constraints <- list(soc,
                      x[1] + x[2] + 3*x[3] >= 1.0,
                      y <= 5)
  obj <- Minimize(3*x[1] + 2*x[2] + x[3])
  expect_x <- matrix(c(-3.874621860638774, -2.129788233677883, 2.33480343377204))
  expect_x <- round(expect_x, 5)
  expect_y <- 5
  var_pairs <- list(list(x, expect_x), list(y, expect_y))
  expect_soc <- list(matrix(2.86560262), matrix(c(2.22062583, 1.22062583, -1.33812252)))
  expect_ineq1 <- 0.7793969212001993
  expect_ineq2 <- 2.865602615049077
  con_pairs <- list(list(constraints[[1]], expect_soc), 
                    list(constraints[[2]], expect_ineq1),
                    list(constraints[[3]], expect_ineq2))
  obj_pair <- list(obj, -13.548638904065102)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

socp_2 <- function() {
  # An (unnecessarily) SOCP-based reformulation of LP_1.
  x <- Variable(nrow = 2, name = "x")
  objective <- Minimize(-4*x[1] - 5*x[2])
  expr <- reshape_expr(x[1] + 2*x[2], c(1, 1))
  constraints <- list(2*x[1] + x[2] <= 3,
                      SOC(Constant(3), expr),
                      x[1] >= 0,
                      x[2] >= 0)
  con_pairs <- list(list(constraints[[1]], 1),
                    list(constraints[[2]], list(matrix(2), matrix(-2))),
                    list(constraints[[3]], 0),
                    list(constraints[[4]], 0))
  var_pairs <- list(list(x, matrix(c(1,1))))
  obj_pair <- list(objective, -9)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

socp_3 <- function(axis) {
  x <- Variable(nrow = 2)
  c <- matrix(c(-1, 2))
  root2 <- sqrt(2)
  u <- rbind(c(1/root2, -1/root2), c(1/root2, 1/root2))
  mat1 <- diag(c(root2, 1/root2)) %*% t(u)
  mat2 <- diag(c(1,1))
  mat3 <- diag(c(0.2, 1.8))
  
  X <- rbind(mat1 %*% x, mat2 %*% x, mat3 %*% x)   # Stack these as rows.
  t <- Constant(matrix(c(1,1,1)))
  objective <- Minimize(t(c) %*% x)
  if(axis == 1) {
    con <- SOC(t, X, axis = 1)
    con_expect <- list(matrix(c(0, 1.16454469e+00, 7.67560451e-01)),
                       rbind(c(0, 0),
                             c(-9.74311819e-01, 6.37872081e-01),
                             c(-1.28440860e-01, 7.56737724e-01)))
  } else {
    con <- SOC(t, t(X), axis = 2)
    con_expect <- list(matrix(c(0, 1.16454469e+00,  7.67560451e-01)),
                       rbind(c(0, -9.74311819e-01, -1.28440860e-01),
                             c(0,  6.37872081e-01,  7.56737724e-01)))
  }
  
  obj_pair <- list(objective, -1.932105)
  con_pairs <- list(list(con, con_expect))
  var_pairs <- list(list(x, matrix(c(0.83666003, -0.54772256))))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

sdp_1 <- function(objective_sense) {
  # Solve "Example 8.3" from Convex Optimization by Boyd & Vandenberghe.
  # Verify (1) optimal objective values, (2) that the dual variable to the PSD constraint
  #  belongs to the correct cone (i.e. the dual variable is itself PSD), and (3) that
  #  complementary slackness holds with the PSD primal variable and its dual variable.
  
  rho <- Variable(4, 4, symmetric = TRUE)
  constraints <- list( 0.6 <= rho[1,2], rho[1,2] <=  0.9,
                       0.8 <= rho[1,3], rho[1,3] <=  0.9,
                       0.5 <= rho[2,4], rho[2,4] <=  0.7,
                      -0.8 <= rho[3,4], rho[3,4] <= -0.4,
                      rho[1,1] == 1, rho[2,2] == 1, rho[3,3] == 1, rho[4,4] == 1, rho %>>% 0)
  
  if(objective_sense == "min") {
    obj <- Minimize(rho[1,4])
    obj_pair <- list(obj, -0.39)
  } else if(objective_sense == "max") {
    obj <- Maximize(rho[1,4])
    obj_pair <- list(obj, 0.23)
  } else
    stop("Unknown objective_sense")
  
  con_pairs <- lapply(constraints, function(c) { list(c, NULL) })
  var_pairs <- list(list(rho, NULL))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

sdp_2 <- function() {
  # Example SDO2 from MOSEK 9.2 documentation.
  X1 <- Variable(2, 2, symmetric = TRUE)
  X2 <- Variable(4, 4, symmetric = TRUE)
  C1 <- rbind(c(1,0), c(0,6))
  A1 <- rbind(c(1,1), c(1,2))
  C2 <- rbind(c(1, -3, 0, 0), c(-3, 2, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 0))
  A2 <- rbind(c(0, 1, 0, 0), c(1, -1, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, -3))
  b <- 23
  k <- -3
  var_pairs <- list(
    list(X1, rbind(c(21.04711571, 4.07709873),
                   c( 4.07709873, 0.7897868))),
    list(X2, rbind(c(5.05366214, -3., 0.,  0.),
                   c(-3., 1.78088676, 0.,  0.),
                   c( 0., 0.,         0.,  0.),
                   c( 0., 0.,         0., -0.)))
  )
  con_pairs <- list(
    list(matrix_trace(A1 %% X1) + matrix_trace(A2 %*% X2) == b, -0.83772234),
    list(X2[1,2] <= k, 11.04455278),
    list(X1 %>>% 0, rbind(c(21.04711571, 4.07709873),
                          c( 4.07709873, 0.7897868))),
    list(X2 %>>% 0, rbind(c(1., 1.68455405, 0., 0.),
                          c(1.68455405, 2.83772234, 0., 0.),
                          c(0.,         0.,         1., 0.),
                          c(0.,         0.,         0., 2.51316702)))
  )
  obj_expr <- Minimize(matrix_trace(C1 %*% X1) + matrix_trace(C2 %*% X2))
  obj_pair <- list(obj_expr, 52.40127214)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

expcone_1 <- function() {
  # min   3*x[1] + 2*x[2] + x[3]
  # s.t.  0.1 <= x[1] + x[2] + x[3] <= 1
  #       x >= 0
  #       x[1] >= x[2]*exp(x[3]/x[2])
  
  x <- Variable(3, 1)
  cone_con <- ExpCone(x[3], x[2], x[1])
  constraints <- list(sum(x) <= 1.0, sum(x) >= 0.1, x >= 0, cone_con)
  obj <- Minimize(3*x[1] + 2*x[2] + x[3])
  obj_pair <- list(obj, 0.23534820622420757)
  expect_eqp <- list(matrix(-1.35348213), matrix(-0.35348211), matrix(0.64651792))
  con_pairs <- list(list(constraints[[1]], 0),
                    list(constraints[[2]], 2.3534821130067614),
                    list(constraints[[3]], matrix(0, nrow = 3, ncol = 1)),
                    list(constraints[[4]], expect_eqp)
                   )
  expect_x <- matrix(c(0.05462721, 0.02609378, 0.01927901), nrow = 3, ncol = 1)
  var_pairs <- list(list(x, expect_x))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

expcone_socp_1 <- function() {
  # A random risk-parity portfolio optimization problem.
  sigma <- rbind(c(1.83, 1.79, 3.22), 
                 c(1.79, 2.18, 3.18), 
                 c(3.22, 3.18, 8.69))
  L <- t(chol(sigma))
  c <- 0.75
  t <- Variable(name = "t")
  x <- Variable(3, name = "x")
  s <- Variable(3, name = "s")
  e <- Constant(matrix(1, nrow = 3))
  
  objective <- Minimize(t - c * t(e) %*% s)
  con1 <- (norm(t(L) %*% x, "2") <= t)
  con2 <- ExpCone(s, e, x)
  
  # SolverTestHelper data.
  obj_pair <- list(objective, 4.0751197)
  var_pairs <- list(
    list(x, matrix(c(0.576079, 0.54315, 0.28037))),
    list(s, matrix(c(-0.55150, -0.61036, -1.27161)))
  )
  con_pairs <- list(
    list(con1, 1.0),
    list(con2, list(matrix(c(-0.75, -0.75, -0.75)),
                    matrix(c(-1.16363, -1.20777, -1.70371)),
                    matrix(c(1.30190, 1.38082, 2.67496)))
        )
    )
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

pcp_1 <- function() {
  # Use a 3D power cone formulation for
  # min 3*x[1] + 2*x[2] + x[3]
  # s.t. norm(x,2) <= y
  #      x[1] + x[2] + 3*x[3] >= 1.0
  #      y <= 5
  
  x <- Variable(3)
  y_square <- Variable()
  epis <- Variable(3)
  constraints <- list(PowCone3D(rep(1, 3), epis, x, Constant(c(0.5, 0.5, 0.5))),
                      sum(epis) <= y_square,
                      x[1] + x[2] + 3*x[3] >= 1.0,
                      y_square <- 25)
  
  obj <- Minimize(3*x[1] + 2*x[2] + x[3])
  expect_x <- matrix(c(-3.874621860638774, -2.129788233677883, 2.33480343377204))
  expect_epis <- expect_x^2
  expect_x <- round(expect_x, 5)
  expect_epis <- round(expect_epis, 5)
  expect_y_square <- 25
  var_pairs <- list(list(x, expect_x),
                    list(epis, expect_epis),
                    list(y_square, expect_y_square))
  expect_ineq1 <- 0.7793969212001993
  expect_ineq2 <- 2.865602615049077/10
  expect_pc <- list(matrix(c(4.30209047, 1.29985494, 1.56211543)),
                    matrix(c(0.28655796, 0.28655796, 0.28655796)),
                    matrix(c(2.22062898, 1.22062899, -1.33811302))
                   )
  con_pairs <- list(list(constraints[[1]], expect_pc),
                    list(constraints[[2]], expect_ineq2),
                    list(constraints[[3]], expect_ineq1),
                    list(constraints[[4]], expect_ineq2))
  obj_pair <- list(obj, -13.548638904065102)
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

pcp_2 <- function() {
  # Reformulate
  #   max  (x**0.2)*(y**0.8) + z**0.4 - x
  #   s.t. x + y + z/2 == 2
  #        x, y, z >= 0
  # Into
  #   max  x3 + x4 - x0
  #   s.t. x0 + x1 + x2 / 2 == 2,
  #        (x0, x1, x3) in Pow3D(0.2)
  #        (x2, 1.0, x4) in Pow3D(0.4)
  
  x <- Variable(3)
  hypos <- Variable(2)
  objective <- Minimize(-sum(hypos) + x[1])
  arg1 <- cbind(x[1], x[3])
  arg2 <- cbind(x[2], 1.0)
  pc_con <- PowCone3D(arg1, arg2, hypos, c(0.2, 0.4))
  expect_pc_con <- list(matrix(1.48466366, 0.24233184),
                        matrix(0.48466367, 0.83801333),
                        matrix(-1, -1))
  con_pairs <- list(
    list(x[1] + x[2] + 0.5*x[3] == 2, 0.4846636697795672),
    list(pc_con, expect_pc_con)
  )
  obj_pair <- list(objective, -1.8073406786220672)
  var_pairs <- list(
    list(x, matrix(c(0.06393515, 0.78320961, 2.30571048))),
    list(hypos, NULL)
  )
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

pcp_3 <- function() {
  w <- Variable(2, 1)
  D <- rbind(c(-1.0856306,   0.99734545),
             c( 0.2829785,  -1.50629471),
             c(-0.57860025,  1.65143654),
             c(-2.42667924, -0.42891263),
             c( 1.26593626, -0.8667404),
             c(-0.67888615, -0.09470897),
             c(1.49138963, -0.638902))   # T-by-N.
  
  # Minimize ||D @ w||_p s.t. 0 <= w, sum(w) == 1.
  # Refer to https://docs.mosek.com/modeling-cookbook/powo.html#p-norm-cones
  p <- 1/0.4
  Tnum <- nrow(D)
  t <- Variable()
  d <- Variable(Tnum, 1)
  ones <- matrix(1, nrow = Tnum, ncol = 1)
  
  powcone <- PowCone3D(d, t*ones, D %*% w, 1/p)
  constraints <- list(sum(w) == 1, w >= 0, powcone, sum(d) == t)
  con_pairs <- list(
    list(constraints[[1]], -1.51430),
    list(constraints[[2]], matrix(c(0,0))),
    list(constraints[[3]], list(
      matrix(c(0.40000935, 0.40000935, 0.40000935, 0.40000935, 0.40000935, 0.40000935, 0.40000935), ncol = 1),
      matrix(c(2.84369172e-03, 1.22657446e-01, 1.12146997e-01, 3.45802205e-01, 2.76327461e-05, 1.27539057e-02, 3.75878155e-03), ncol = 1),
      matrix(c(-0.04031276, 0.38577107, -0.36558292, 0.71847219, 0.00249992, 0.09919715, -0.04765863))
    )),
    list(constraints[[4]], 0.40000935)
  )
  
  univar_obj <- function(w0) {
    v <- D[,1]*w0 + D[,2]*(1 - w0)
    return(sum(v^p)^(1/p))
  }
  
  univar_res <- optim(matrix(0.4), univar_obj, lower = 0, upper = 1, abstol = 1e-16)   # Is abstol here equivalent to tol in scipy.optimize.minimize?
  w_opt <- matrix(c(res$par, 1 - res$par))
  
  univar_res_fun <- univar_obj(univar_res$x)
  obj_pair <- list(Minimize(t), univar_res_fun)
  var_pairs <- list(list(d, matrix(c(7.17144981e-03, 3.09557056e-01, 2.83038570e-01, 8.72785905e-01, 6.92995408e-05, 3.21904516e-02, 9.48918352e-03))),
                    list(w, w_opt),
                    list(t, matrix(univar_res_fun)))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

mi_lp_0 <- function() {
  x <- Variable(2)
  bool_var <- Variable(boolean = TRUE)
  con_pairs <- list(list(x == bool_var, NULL),
                    list(bool_var == 0, NULL))
  obj_pair <- list(Minimize(norm(x, 1) + 1.0), 1)
  var_pairs <- list(list(x, matrix(c(0,0))),
                    list(bool_var, 0))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

mi_lp_1 <- function() {
  x <- Variable(2, name = "x")
  boolvar <- Variable(boolean = TRUE)
  intvar <- Variable(integer = TRUE)
  objective <- Minimize(-4*x[1] - 5*x[2])
  constraints <- list(2*x[1] + x[2] <= intvar,
                      x[1] + 2*x[2] <= 3*boolvar,
                      x >= 0,
                      intvar == 3*boolvar,
                      intvar == 3)
  obj_pair <- list(objective, -9)
  var_pairs <- list(list(x, matrix(c(1,1))),
                    list(boolvar, 1),
                    list(intvar, 3))
  con_pairs <- lapply(constraints, function(c) { list(c, NULL) })
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}


