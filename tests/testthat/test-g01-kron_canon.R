context("test-g01-kron_canon")
TOL <- 1e-6

# The Kronecker product of matrices M, N is:
#     kron(M, N) = [M[0,0] * N   , ..., M[0, end] * N  ]
#                  [M[1,0] * N   , ..., M[1, end] * N  ]
#                  ...
#                  [M[end, 0] * N, ..., M[end, end] * N]

make_kron_prob <- function(z_dims, c_dims, param, var_left, seed) {
  # Construct random nonnegative matrices (C, L) of shapes
  # (c_dims, z_dims) respectively. Define an optimization
  # problem with a matrix variable of shape z_dims:
  #   
  #   min sum(Z)
  #   s.t.  kron(Z, C) >= kron(L, C)   ---   if var_left is True
  #         kron(C, Z) >= kron(C, L)   ---   if var_left is False
  #         Z >= 0
  # 
  # Regardless of whether var_left is True or False, the optimal
  # solution to that problem is Z = L.
  # 
  # If param is TRUE, then C is defined as a CVXR Parameter.
  # If param is FALSE, then C is a CVXR Constant.
  # 
  # A small remark: the constraint that Z >= 0 is redundant.
  # It's there because it's easier to set break points that distinguish
  # objective canonicalization and constraint canonicalization
  # when there's more than one constraint.
  
  set.seed(seed)
  
  C_value <- matrix(runif(prod(c_dims)), nrow = c_dims[1], ncol = c_dims[2])
  C_value <- round(C_value, digits = 2)
  if(param) {
    C <- new("Parameter", dim = c_dims)
    value(C) <- C_value
  } else
    C <- Constant(C_value)
  
  Z <- new("Variable", dim = z_dims)
  L <- matrix(runif(prod(z_dims)), nrow = z_dims[1], ncol = z_dims[2])
  L <- round(L, digits = 2)
  
  if(var_left) {
    constraints <- list(kronecker(Z, C) >= kronecker(L, C), Z >= 0)
    # The cvxcore function get_kronl_mat doesn't work when C is a Parameter.
    # We get around this by having kron be non-dpp, but this comes at
    # the price of eliminating the speed benefit of using Parameter objects.
    # We'll eventually need to extend get_kronl_mat so that it supports
    # Parameters. Until then, I'll make a note that tests here DO PASS
    # with the existing get_kronl_mat implementation if we use the following
    # constraints: list(kronecker(Z - L, C) >= 0, Z >= 0).
  } else
    constraints <- list(kronecker(C, Z) >= kronecker(C, L), Z >= 0)
  
  obj_expr <- sum(Z)
  prob <- Problem(Minimize(obj_expr), constraints)
  return(list(Z, C, L, prob))
}

TestKronRightVar.C_DIMS <- list(c(1,1), c(2,1), c(1,2), c(2,2))
TestKronLeftVar.C_DIMS <- list(c(1,1), c(2,1), c(1,2), c(2,2))

symvar_kronl <- function(param) {
  # Use a symmetric matrix variable.
  X <- Variable(2, 2, symmetric = TRUE)
  b_val <- matrix(1.5)
  if(param) {
    b <- Parameter(1,1)
    value(b) <- b_val
  } else
    b <- Constant(b_val)
  
  L <- rbind(c(0.5, 1), c(2, 3))
  U <- rbind(c(10, 11), c(12, 13))
  kronX <- kronecker(X, b)   # Should be equal to X.
  
  objective <- Minimize(sum(flatten(X)))
  constraints <- list(U >= kronX, kronX >= L)
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  
  expect_equal(result$getValue(X), rbind(c(0.5, 2), c(2, 3)) / 1.5)
  objective <- Maximize(sum(flatten(X)))
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  expect_equal(result$getValue(X), rbind(c(10, 11), c(11, 13)) / 1.5)
}

scalar_kronl <- function(param) {
  y <- Variable(1,1)
  A_val <- rbind(c(1,2), c(3,4))
  L <- rbind(c(0.5,1), c(2,3))
  U <- rbind(c(10,11), c(12,13))
  if(param) {
    A <- Parameter(2,2)
    value(A) <- A_val
  } else
    A <- Constant(A_val)
  
  krony <- kronecker(y, A)   # Should be equal to y*A.
  constraints <- list(U >= krony, krony >= L)
  
  objective <- Minimize(y)
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  expect_equal(result$getValue(y), matrix(max(L / A_val)), tolerance = TOL)
  
  objective <- Maximize(y)
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  expect_equal(result$getValue(y), matrix(min(U / A_val)), tolerance = TOL)
}

test_that("TestKronRightVar: test gen kronr param", {
  z_dims <- c(2,2)
  for(c_dims in TestKronRightVar.C_DIMS) {
    tmp <- make_kron_prob(z_dims, c_dims, param = TRUE, var_left = FALSE, seed = 0)
    Z <- tmp[[1]]
    C <- tmp[[2]]
    L <- tmp[[3]]
    prob <- tmp[[4]]
    
    result <- solve(prob, solver = "ECOS", abstol = 1e-8, reltol = 1e-8)
    expect_equal(result$status, OPTIMAL)
    con_viols <- result$getViolation(prob@constraints[[1]])
    expect_lte(max(con_viols), 1e-4)
    expect_equal(result$getValue(Z), L, tolerance = 1e-4)
  }
})

test_that("TestKronRightVar: test gen kronr const", {
  z_dims <- c(2,2)
  for(c_dims in TestKronRightVar.C_DIMS) {
    tmp <- make_kron_prob(z_dims, c_dims, param = FALSE, var_left = FALSE, seed = 0)
    Z <- tmp[[1]]
    C <- tmp[[2]]
    L <- tmp[[3]]
    prob <- tmp[[4]]
    
    result <- solve(prob, solver = "ECOS", abstol = 1e-8, reltol = 1e-8)
    expect_equal(result$status, OPTIMAL)
    con_viols <- result$getViolation(prob@constraints[[1]])
    expect_lte(max(con_viols), 1e-4)
    expect_equal(result$getValue(Z), L, tolerance = 1e-4)
  }
})

test_that("TestKronLeftVar: test symvar kronl param", {
  symvar_kronl(param = TRUE)
})

test_that("TestKronLeftVar: test symvar kronl const", {
  symvar_kronl(param = FALSE)
})

test_that("TestKronLeftVar: test scalar kronl param", {
  scalar_kronl(param = TRUE)
})

test_that("TestKronLeftVar: test scalar kronl const", {
  scalar_kronl(param = FALSE)
})

test_that("test gen kronl param", {
  z_dims <- c(2,2)
  for(c_dims in TestKronLeftVar.C_DIMS) {
    tmp <- make_kron_prob(z_dims, c_dims, param = TRUE, var_left = TRUE, seed = 0)
    Z <- tmp[[1]]
    C <- tmp[[2]]
    L <- tmp[[3]]
    prob <- tmp[[4]]
    
    result <- solve(prob, solver = "ECOS", abstol = 1e-8, reltol = 1e-8)
    expect_equal(result$status, OPTIMAL)
    con_viols <- result$getViolation(prob@constraints[[1]])
    expect_lte(max(con_viols), 1e-4)
    expect_equal(result$getValue(Z), L, tolerance = 1e-4)
  }
})

test_that("test gen kronl const", {
  z_dims <- c(2,2)
  for(c_dims in TestKronLeftVar.C_DIMS) {
    tmp <- make_kron_prob(z_dims, c_dims, param = FALSE, var_left = TRUE, seed = 0)
    Z <- tmp[[1]]
    C <- tmp[[2]]
    L <- tmp[[3]]
    prob <- tmp[[4]]
    
    result <- solve(prob, solver = "ECOS", abstol = 1e-8, reltol = 1e-8)
    expect_equal(result$status, OPTIMAL)
    con_viols <- result$getViolation(prob@constraints[[1]])
    expect_lte(max(con_viols), 1e-4)
    expect_equal(result$getValue(Z), L, tolerance = 1e-4)
  }
})
