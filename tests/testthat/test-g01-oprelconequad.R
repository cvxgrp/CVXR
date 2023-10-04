context("test-g01-oprelconequad")
TOL <- 1e-5

n <- 3

a <- Variable(n, pos = TRUE)
b <- Variable(n, pos = TRUE)

a_lower <- cumsum(runif(n))
a_upper <- a_lower + 0.05*runif(n)
b_lower <- cumsum(runif(n))
b_upper <- b_lower + 0.05*runif(n)

base_cons <- list(a_lower <= a, a <= a_upper,
                  b_lower <= b, b <= b_upper)

installed_solvers_vec <- installed_solvers()

TestOpRelConeQuad.get_solver <- function() {
  if("MOSEK" %in% installed_solvers_vec)
    return("MOSEK")
  else if("CVXOPT" %in% installed_solvers_vec)
    return("CVXOPT")
  else if("COPT" %in% installed_solvers_vec)
    return("COPT")
  else
    stop("No viable solver installed")
}

solver <- TestOpRelConeQuad.get_solver()

TestOpRelConeQuad.Dop_commute <- function(ac, bc, U) {
  D <- diag(ac*log(ac/bc))
  if(is.complex(U))
    out <- U %*% D %*% t(Conj(U))
  else
    out <- U %*% D %*% t(U)
  return(out)
}

TestOpRelConeQuad.sum_rel_entr_approx <- function(ae, be, apx_m, apx_k) {
  nc <- size(ae)
  if(nc != size(be))
    stop("ae and be must be the same size")
  epi_vec <- Variable(nc)
  con <- RelEntrConeQuad(ae, be, epi_vec, apx_m, apx_k)
  objective <- Minimize(sum(epi_vec))
  return(list(objective, con))
}

oprelcone_1 <- function(apx_m, apx_k, real) {
  # These tests construct two matrices that commute (imposing all eigenvectors equal)
  # and then use the fact that: T=Dop(A, B) for (A, B, T) in OpRelEntrConeQuad
  # i.e. T >> Dop(A, B) for an objective that is an increasing function of the
  # eigenvalues (which we here take to be the trace), we compute the reference
  # objective value as tr(Dop) whose correctness can be seen by writing out
  # tr(T)=tr(T-Dop)+tr(Dop), where tr(T-Dop)>=0 because of PSD-ness of (T-Dop),
  # and at optimality we have (T-Dop)=0 (the zero matrix of corresponding size)
  # For the case that the input matrices commute, Dop takes on a particularly
  # simplified form, i.e.: U @ diag(a * log(a/b)) @ U^{-1} (which is implemented
  # in the Dop_commute method above).
  
  # Compute the expected optimal solution.
  tmp <- TestOpRelConeQuad.sum_rel_entr_approx(a, b, apx_m, apx_k)
  temp_obj <- tmp[[1]]
  temp_con <- tmp[[2]]
  temp_constraints <- base_cons
  temp_constraints <- c(temp_constraints, temp_con)
  temp_prob <- Problem(temp_obj, temp_constraints)
  result <- solve(temp_prob)
  expect_a <- result$getValue(a)
  expect_b <- result$getValue(b)
  expect_objective <- result$getValue(temp_obj)
  
  # Next: Create a matrix representation of the same problem, using operator
  # relative entropy.
  nc <- n
  if(real) {
    randmat <- matrix(rnorm(nc^2), nrow = nc, ncol = nc)
    QR_res <- qr(randmat)
    U <- qr.Q(QR_res)
    A <- symmetric_wrap(U %*% diag(a) %*% t(U))
    B <- symmetric_wrap(U %*% diag(b) %*% t(U))
    Tvar <- Variable(n, n, symmetric = TRUE)
  } else {
    randmat <- 1i*matrix(rnorm(nc^2), nrow = nc, ncol = nc)
    randmat <- randmat + matrix(rnorm(nc^2), nrow = nc, ncol = nc)
    QR_res <- qr(randmat)
    U <- qr.Q(QR_res)
    A <- hermitian_wrap(U %*% diag(a) %*% t(Conj(U)))
    B <- hermitian_wrap(U %*% diag(b) %*% t(Conj(U)))
    Tvar <- Variable(n, n, hermitian = TRUE)
  }
  
  main_con <- OpRelEntrConeQuad(A, B, Tvar, apx_m, apx_k)
  obj <- Minimize(matrix_trace(Tvar))
  expect_T <- TestOpRelConeQuad.Dop_commute(expect_a, expect_b, U)
  
  # Define the SolverTestHelper object.
  con_pairs <- lapply(base_cons, function(con) { list(con, NULL) })
  con_pairs <- c(con_pairs, list(list(main_con, NULL)))
  obj_pairs <- list(obj, expect_objective)
  var_pairs <- list(list(Tvar, expect_T))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("test oprelcone 1, m = 1, k = 3, real", {
  sth <- oprelcone_1(1, 3, TRUE)
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_objective(sth, result, tolerance = 1e-3)
})

test_that("test oprelcone 1, m = 3, k = 1, real", {
  sth <- oprelcone_1(3, 1, TRUE)
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_objective(sth, result, tolerance = 1e-3)
})

test_that("test oprelcone 1, m = 4, k = 4, real", {
  sth <- oprelcone_1(4, 4, TRUE)
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_objective(sth, result, tolerance = 1e-3)
})

test_that("test oprelcone 1, m = 1, k = 3, complex", {
  sth <- oprelcone_1(1, 3, FALSE)
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_objective(sth, result, tolerance = 1e-3)
})

test_that("test oprelcone 1, m = 3, k = 1, complex", {
  sth <- oprelcone_1(3, 1, FALSE)
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-3)
  verify_objective(sth, result, tolerance = 1e-3)
})

oprelcone_2 <- function() {
  # This test uses the same idea from the tests with commutative matrices,
  # instead, here, we make the input matrices to Dop, non-commutative,
  # the same condition as before i.e. T=Dop(A, B) for (A, B, T) in OpRelEntrConeQuad
  # (for an objective that is an increasing function of the eigenvalues) holds,
  # the difference here then, is in how we compute the reference values, which
  # has been done by assuming correctness of the original CVXQUAD matlab implementation
  
  nc <- 4
  m <- 3
  k <- 3
  
  # Generate two sets of linearly orthogonal vectors.
  # Each to be set as the eigenvectors of a particular input matrix to Dop.
  U1 <- rbind(c(-0.05878522, -0.78378355, -0.49418311, -0.37149791),
              c( 0.67696027, -0.25733435,  0.59263364, -0.35254672),
              c( 0.43478177,  0.53648704, -0.54593428, -0.47444939),
              c( 0.59096015, -0.17788771, -0.32638042,  0.71595942))
  U2 <- rbind(c(-0.42499169,  0.6887562,   0.55846178, 0.18198188),
              c(-0.55478633, -0.7091174,   0.3884544,  0.19613213),
              c(-0.55591804,  0.14358541, -0.72444644, 0.38146522),
              c( 0.4500548,  -0.04637494,  0.11135968, 0.88481584))
  
  a_diag <- Variable(n, pos = TRUE)
  b_diag <- Variable(n, pos = TRUE)
  A <- U1 %*% diag(a_diag) %*% t(U1)
  B <- U2 %*% diag(b_diag) %*% t(U2)
  Tvar <- Variable(n, n, symmetric = TRUE)
  
  ac_lower <- c(0.40683013, 1.34514597, 1.60057343, 2.13373667)
  ac_upper <- c(1.36158501, 1.61289351, 1.85065805, 3.06140939)
  bc_lower <- c(0.06858235, 0.36798274, 0.95956627, 1.16286541)
  bc_upper <- c(0.70446555, 1.16635299, 1.46126732, 1.81367755)
  
  con1 <- OpRelEntrConeQuad(A, B, Tvar, m, k)
  con2 <- ac_lower <= a_diag
  con3 <- a_diag <= ac_upper
  con4 <- bc_lower <= b_diag
  con5 <- b_diag <= bc_upper
  con_pairs <- list(list(con1, NULL),
                    list(con2, NULL),
                    list(con3, NULL),
                    list(con4, NULL),
                    list(con5, NULL))
  obj <- Minimize(matrix_trace(Tvar))
  
  expect_obj <- 1.85476
  expect_T <- rbind(c( 0.49316819, 0.20845265, 0.60474713, -0.5820242),
                    c( 0.20845265, 0.31084053, 0.2264112,  -0.8442255),
                    c( 0.60474713, 0.2264112,  0.4687153,  -0.85667283),
                    c(-0.5820242, -0.8442255, -0.85667283,  0.58206723))
  
  obj_pair <- list(obj, expect_obj)
  var_pairs <- list(list(Tvar, expect_T))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

test_that("test oprelcone 2", {
  sth <- oprelcone_2()
  result <- solve(sth, solver)
  verify_primal_values(sth, result, tolerance = 1e-2)
  verify_objective(sth, result, tolerance = 1e-2)
})











