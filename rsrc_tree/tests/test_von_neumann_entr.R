context("test_von_neumann_entr")

get_solver_args <- function() {
  if("MOSEK" %in% installed_solvers())
    SOLVE_ARGS <- list(solver = "MOSEK", verbose =  TRUE)
  else
    SOLVE_ARGS <- list(solver = "SCS", eps = 1e-6, max_iters = 5e5, verbose = TRUE)
  return(SOLVE_ARGS)
}

SOLVE_ARGS <- get_solve_args()

make_test_1 <- function(complex) {
  # Enforce an upper bound of 0.8 on trace(N);
  # Expect N's unspecified eigenvalue to be 0.2.
  n <- 3
  
  if(complex) {
    N <- Variable(n, n, hermitian = TRUE)
    V12 <- matrix(rnorm(n*n), nrow = n, ncol = n) + 1i*matrix(rnorm(n*n), nrow = n, ncol = n)
  } else {
    N <- Variable(n, n, PSD = TRUE)
    V12 <- matrix(rnorm(n*n), nrow = n, ncol = n)
  }
  
  V12 <- qr.Q(qr(V12))[,1:2]
  mu12 <- matrix(c(0.5,0.1))
  trace_bound <- 0.2
  cons1 <- N %*% V12 == V12 * mu12
  cons2 <- matrix_trace(N) <= trace_bound
  objective <- Maximize(von_neumann_entr(N))
  
  V3 <- matrix(onb_for_orthogonal_complement(V12), nrow = n, ncol = 1, byrow = TRUE)
  mu2 <- trace_bound - sum(mu12)
  expect_mu <- c(mu12, mu3)
  expect_V <- cbind(V12, V3)
  if(complex)
    expect_N <- (expect_V*expect_mu) %*% t(Conj(expect_V))
  else
    expect_N <- (expect_V*expect_mu) %*% t(expect_V)
  expect_obj <- value(sum(entr(expect_mu)))
  
  obj_pair <- list(objective, expect_obj)
  con_pairs <- list(list(cons1, NA_real_), list(cons2, NA_real_))
  var_pairs <- list(list(N, expect_N))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

make_test_2 <- function(quad_approx) {
  # Enforce a lower bound of 0.9 on trace(N);
  # Expect N's unspecified eigenvalue to be 0.4.
  n <- 3
  N <- Variable(n, n, PSD = TRUE)
  V12 <- rbind(c(-0.12309149, 0.90453403),
               c(-0.49236596, 0.30151134),
               c(-0.86164044, -0.30151134))
  mu12 <- matrix(c(0.3, 0.2))
  trMin <- 0.9
  cons1 <- N %*% V12 == V12*mu12
  cons2 <- matrix_trace(N) >= trMin
  if(quad_approx)
    objective <- Maximize(von_neumann_entr(N, c(5,5)))
  else
    objective <- Maximize(von_neumann_entr(N))
  
  V3 <- matrix(onb_for_orthogonal_complement(V12), nrow = n, ncol = 1, byrow = TRUE)
  mu3 <- trMin - sum(mu12)
  expect_mu <- c(mu12, mu3)
  expect_V <- cbind(V12, V3)
  expect_N <- (expect_V*expect_mu) %*% t(expect_V)
  expect_obj <- value(sum(entr(expect_mu)))
  
  obj_pair <- list(objective, expect_obj)
  con_pairs <- list(list(cons1, NA_real_), list(cons2, NA_real_))
  var_pairs <- list(list(N, expect_N))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

sum_entr_approx <- function(a, apx_m, apx_k) {
  n <- size(a)
  epi_vec <- Variable(n)
  b <- Constant(rep(1,n))
  con <- RelEntrConeQuad(a, b, epi_vec, apx_m, apx_k)
  objective <- Minimize(sum(epi_vec))
  return(list(objective, con))
}

make_test_3 <- function(quad_approx = FALSE, real = FALSE) {
  set.seed(0)
  
  ###################################################
  #
  #   Construct matrix/vector coefficient data
  #
  ###################################################
  
  apx_m <- 2
  apx_k <- 2
  A1_real <- rbind(c(8.38972, 1.02671, 0.87991),
                   c(1.02671, 8.41455, 7.31307),
                   c(0.87991, 7.31307, 2.35915))
  A2_real <- rbind(c(6.92907, 4.37713, 5.11915),
                   c(4.37713, 7.96725, 4.42217),
                   c(5.11915, 4.42217, 2.72919))
  
  if(real) {
    U <- diag(3)
    A1 <- A1_real
    A2 <- A2_real
  } else {
    randmat <- 1i*matrix(rnorm(9), nrow = 3, ncol = 3)
    randmat <- randmat + matrix(rnorm(9), nrow = 3, ncol = 3)
    U <- qr.Q(qr(randmat))
    A1 <- U %*% A1_real %*% t(Conj(U))
    A2 <- U %*% A2_real %*% t(Conj(U))
  }
  b <- matrix(c(19.16342, 17.62551))
  
  ###################################################
  #
  #   define and solve a reference problem
  #
  ###################################################
  diag_X <- Variable(3)
  ref_X <- diag(diag_X)
  if(real) {
    ref_cons <- list(matrix_trace(A1 %*% ref_X) == b[1],
                     matrix_trace(A2 %*% ref_X) == b[2])
  } else {
    conjugated_X <- U %*% ref_X %*% t(Conj(U))
    ref_cons <- list(matrix_trace(A1 %*% conjugated_X) == b[1],
                     matrix_trace(A2 %*% conjugated_X) == b[2])
  }
  
  if(quad_approx) {
    tmp <- sum_entr_approx(diag_X, apx_m, apx_k)
    ref_objective <- tmp[[1]]
    con <- tmp[[2]]
    ref_cons <- c(ref_cons, con)
  } else
    ref_objective <- Minimize(-sum(entr(diag_X)))
  ref_prob <- Problem(ref_objective, ref_cons)
  ref_obj_result <- solve(ref_prob)
  ref_obj_val <- ref_obj_result$value
  
  ###################################################
  #
  #   define a new problem that is equivalent to the
  #   reference, but makes use of von_neumann_entr.
  #
  ###################################################
  if(real) {
    N <- Variable(3, 3, PSD = TRUE)
    cons <- list(matrix_trace(A1 %*% N) == b[1],
                 matrix_trace(A2 %*% N) == b[2],
                 N - diag(diag(N)) == 0)
    expect_N <- value(ref_X)
  } else {
    N <- Variable(3, 3, hermitian = TRUE)
    aconj_N <- t(Conj(U)) %*% N %*% U
    cons <- list(matrix_trace(A1 %*% N) == b[1],
                 matrix_trace(A2 %*% N) == b[2],
                 aconj_N - diag(diag(aconj_N)) == 0)
    expect_N <- value(conjugated_X)
  }
  
  if(quad_approx)
    objective <- Minimize(-von_neumann_entr(N, c(apx_m, apx_k)))
  else
    objective <- Minimize(-von_neumann_entr(N))
  
  ###################################################
  #
  #   construct and return the SolverTestHelper
  #
  ###################################################
  obj_pair <- list(objective, ref_obj_val)
  var_pairs <- list(list(N, expect_N))
  con_pairs <- lapply(cons, function(con) { list(con, NA_real_) })
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("test 1 real", {
  sth <- make_test_1(FALSE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 1 complex", {
  sth <- make_test_1(TRUE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 2 exact", {
  sth <- make_test_2(FALSE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 2 approx", {
  sth <- make_test_2(TRUE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 3 exact real", {
  sth <- make_test_3(quad_approx = FALSE, real = TRUE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 3 approx real", {
  sth <- make_test_3(quad_approx = TRUE, real = TRUE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_primal_feasibility(sth, result, tolerance = 1e-3)
})

test_that("test 3 exact complex", {
  sth <- make_test_3(quad_approx = FALSE, real = FALSE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
})

test_that("test 3 approx complex", {
  sth <- make_test_3(quad_approx = TRUE, real = FALSE)
  result <- do.call("solve", c(list(object = sth), SOLVE_ARGS))
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
  check_primal_feasibility(sth, result, tolerance = 1e-3)
})
