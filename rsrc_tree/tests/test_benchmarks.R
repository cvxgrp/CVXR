context("test_benchmarks")
TOL <- 1e-6

benchmark <- function(func, iters = 1, name = NA_character_) {
  vals <- c()
  for(i in seq(iters)) {
    start <- Sys.time()
    func()
    vals <- c(vals, Sys.time() - start)
  }
  
  if(is.na(name))
    name <- deparse(substitute(func))
  print(paste(name, ": avg = ", sum(vals)/length(vals), " s, std = ", sd(vals), " s (", iters, " iterations)", sep = ""))
}

test_that("test diffcp sdp example", {
  randn_symm <- function(n) {
    A <- matrix(rnorm(n^2), nrow = n, ncol = n)
    return((A + t(A))/2)
  }
  
  randn_psd <- function(n) {
    A <- 1/10 * matrix(rnorm(n^2), nrow = n, ncol = n)
    return(A %*% t(A))
  }
  
  n <- 100
  p <- 100
  C <- randn_psd(n)
  As <- lapply(1:p, function(i) { randn_symm(n) })
  Bs <- matrix(rnorm(p))
  
  diffcp_sdp <- function() {
    X <- Variable(n, n, PSD = TRUE)
    objective <- matrix_trace(C %*% X)
    constraints <- lapply(1:p, function(i) { As[[i]] %*% X == Bs[[i]] })
    problem <- Problem(Minimize(objective), constraints)
    data <- get_problem_data(problem, "SCS")
  }
  
  benchmark(diffcp_sdp, iters = 1)
})

test_that("test tv inpainting", {
  Uorig <- array(rnorm(512*512*3), dim = c(512, 512, 3))
  Udim <- dim(Uorig)
  rows <- Udim[1]
  cols <- Udim[2]
  colors <- Udim[3]
  known <- array(rep(0, rows*cols*colors), dim = c(rows, cols, colors))
  
  for(i in 1:rows) {
    for(j in 1:cols) {
      if(runif(1) > 0.7) {
        for(k in 1:colors)
          known[i, j, k] <- 1
      }
    }
  }
  
  tv_inpainting <- function() {
    Ucorr <- known*Uorig   # This is elementwise multiplication of R arrays.
    variables <- list()
    constraints <- list()
    for(i in 1:colors) {
      U <- Variable(rows, cols)
      variables <- c(variables, U)
      constraints <- c(constraints, multiply(known[,,i], U) == multiply(known[,,i], Ucorr[,,i]))
    }
    problem <- Problem(Minimize(do.call("tv", variables)), constraints)
    data <- get_problem_data(problem, "SCS")
  }
  
  benchmark(tv_inpainting, iters = 1)
})

test_that("test least squares", {
  m <- 20
  n <- 15
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m))
  
  least_squares <- function() {
    x <- Variable(n)
    cost <- sum_squares(A %*% x - b)
    data <- get_problem_data(Problem(Minimize(cost)), "OSQP")
  }
  
  benchmark(least_squares, iters = 1)
})

test_that("test qp", {
  m <- 15
  n <- 10
  p <- 5
  
  P <- matrix(rnorm(n^2), nrow = n, ncol = n)
  P <- t(P) %*% P
  q <- matrix(rnorm(n))
  G <- matrix(rnorm(m*n), nrow = m, ncol = n)
  h <- G %*% matrix(rnorm(n))
  A <- matrix(rnorm(p*n), nrow = p, ncol = n)
  b <- matrix(rnorm(p))
  
  qp <- function() {
    x <- Variable(n)
    get_problem_data(Problem(Minimize((1/2)*quad_form(x, P) + t(q) %*% x), list(G %*%  x <= h, A %*% x == b)), "OSQP")
  }
  
  benchmark(qp, iters = 1)
})

test_that("test cone matrix stuffing with many constraints", {
  m <- 2000
  n <- 2000
  
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  C <- matrix(runif(floor(m/2)))
  b <- matrix(rnorm(m))
  
  x <- Variable(n)
  cost <- sum(A %*% x)
  
  cons1 <- lapply(1:floor(m/2), function(i) { C[i]*x[i] <= b[i] })
  cons2 <- lapply(1:floor(m/2), function(i) { C[i]*x[floor(m/2) + i] == b[floor(m/2) + i] })
  constraints <- c(cons1, cons2)
  
  problem <- Problem(Minimize(cost), constraints)
  
  cone_matrix_stuffing_with_many_constraints <- function() {
    perform(ConeMatrixStuffing(), problem)
  }
  
  benchmark(cone_matrix_stuffing_with_many_constraints, iters = 1)
})

test_that("test parametrized cone matrix stuffing with many constraints", {
  print("Skipping test. This benchmark takes too long")
  return()
  
  m <- 2000
  n <- 2000
  
  A <- Parameter(m, n)
  C <- Parameter(floor(m/2))
  b <- Parameter(m)
  
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(C) <- matrix(runif(floor(m/2)))
  value(b) <- matrix(rnorm(m))
  
  x <- Variable(n)
  cost <- sum(A %*% x)
  
  cons1 <- lapply(1:floor(m/2), function(i) { C[i]*x[i] <= b[i] })
  cons2 <- lapply(1:floor(m/2), function(i) { C[i]*x[floor(m/2) + i] == b[floor(m/2) + i] })
  constraints <- c(cons1, cons2)
  
  problem <- Problem(Minimize(cost), constraints)
  
  parametrized_cone_matrix_stuffing <- function() {
    perform(ConeMatrixStuffing(), problem)
  }
  
  benchmark(parameterized_cone_matrix_stuffing, iters = 1)
})

test_that("test small cone matrix stuffing", {
  m <- 200
  n <- 200
  
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  C <- matrix(runif(floor(m/2)))
  b <- matrix(rnorm(m))
  
  x <- Variable(n)
  cost <- sum(A %*% x)
  
  cons1 <- lapply(1:floor(m/2), function(i) { C[i]*x[i] <= b[i] })
  cons2 <- lapply(1:floor(m/2), function(i) { C[i]*x[floor(m/2) + i] == b[floor(m/2) + i] })
  constraints <- c(cons1, cons2)
  
  problem <- Problem(Minimize(cost), constraints)
  
  small_cone_matrix_stuffing <- function() {
    perform(ConeMatrixStuffing(), problem)
  }
  
  benchmark(small_cone_matrix_stuffing, iters = 10)
})

test_that("test small parameterized cone matrix stuffing", {
  # Failing in Windows CI - potentially memory leak.
  m <- 200
  n <- 200
  
  A <- Parameter(m, n)
  C <- Parameter(floor(m/2))
  b <- Parameter(m)
  
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(C) <- matrix(runif(floor(m/2)))
  value(b) <- matrix(rnorm(m))
  
  x <- Variable(n)
  cost <- sum(A %*% x)
  
  cons1 <- lapply(1:floor(m/2), function(i) { C[i]*x[i] <= b[i] })
  cons2 <- lapply(1:floor(m/2), function(i) { C[i]*x[floor(m/2) + i] == b[floor(m/2) + i] })
  constraints <- c(cons1, cons2)
  
  problem <- Problem(Minimize(cost), constraints)
  
  small_parametrized_cone_matrix_stuffing <- function() {
    perform(ConeMatrixStuffing(), problem)
  }
  
  benchmark(small_parameterized_cone_matrix_stuffing, iters = 1)
})

test_that("test small lp", {
  m <- 200
  n <- 200
  
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m))
  c <- matrix(runif(n))
  
  x <- Variable(n)
  cost <- t(c) %*% x
  constraints <- list(A %*% x <= b)
  problem <- Problem(Minimize(cost), constraints)
  
  small_lp <- function() {
    get_problem_data(problem, "SCS")
  }
  
  benchmark(small_lp, iters = 1)
  benchmark(small_lp, iters = 1, name = "small_lp_second_time")
})

test_that("test small parameterized lp", {
  # Failing in Windows CI - potentially memory leak.
  m <- 200
  n <- 200
  
  A <- Parameter(m, n)
  b <- Parameter(m)
  c <- Parameter(n)
  
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  value(c) <- matrix(runif(n))
  
  x <- Variable(n)
  cost <- t(c) %*% x
  constraints <- list(A %*% x <= b)
  problem <- Problem(Minimize(cost), constraints)
  
  small_parameterized_lp <- function() {
    get_problem_data(problem, "SCS")
  }
  
  benchmark(small_parameterized_lp, iters = 1)
  benchmark(small_parameterized_lp, iters = 1, name = "small_parameterized_lp_second_time")
})

test_that("test parameterized qp", {
  # Test speed of first solve with QP codepath and SOCP codepath.
  m <- 150
  n <- 100
  
  set.seed(1)
  A <- Parameter(m, n)
  b <- Parameter(m)
  
  x <- Variable(n)
  objective <- Minimize(sum_squares(A %*% x - b))
  constraints <- list(0 <= x, x <= 1)
  prob <- Problem(objective, constraints)
  
  start <- Sys.time()
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  result <- solve(prob, solver = "ECOS")
  end <- Sys.time()
  
  print("Conic canonicalization")
  print(paste("(ECOS) solver time:", result$solver_stats$solve_time))
  print(paste("CVXR time:", (end - start) - result$solver_stats$solve_time))
  
  set.seed(1)
  A <- Parameter(m, n)
  b <- Parameter(m)
  
  x <- Variable(n)
  objective <- Minimize(sum_squres(A %*% x - b))
  constraints <- list(0 <= x, x <= 1)
  prob <- Problem(objective, constraints)
  
  start <- Sys.time()
  value(A) <- matrix(rnorm(m*n), nrow = m, ncol = n)
  value(b) <- matrix(rnorm(m))
  result <- solve(prob, solver = "OSQP")
  end <- Sys.time()
  
  print("Quadratic canonicalization")
  print(paste("(OSQP) solver time:", result$solver_stats$solve_time))
  print(paste("CVXR time:", (end - start) - result$solver_stats$solve_time))
})

test_that("test issue 1668 slow pruning", {
  # Regression test for https://github.com/cvxpy/cvxpy/issues/1668
  # Pruning matrices caused order-of-magnitude slowdowns in compilation times.
  s <- 2000
  t <- 10
  x <- seq(-100.0, 100.0, length.out = s)
  rows <- 50
  var <- Variable(rows, t)
  
  A <- matrix(do.call("rbind", lapply(1:t, function(i) { x })), nrow = t, ncol = nrow(x), byrow = TRUE)
  B <- matrix(rep(x, rows), nrow = rows, ncol = s, byrow = TRUE)
  cost <- sum_squares(var %*% A - B)
  objective <- Minimize(cost)
  problem <- Problem(objective)
  
  start <- Sys.time()
  data <- get_problem_data(problem, "ECOS", verbose = TRUE)
  end <- Sys.time()
  
  print("Issue #1668 regression test")
  print(paste("Compilation test:", end - start))
})







