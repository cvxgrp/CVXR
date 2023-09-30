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

