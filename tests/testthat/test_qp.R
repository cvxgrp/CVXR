context("test-g01-qp")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")
w <- Variable(5, name = "w")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

slope <- Variable(1, name = "slope")
offset <- Variable(1, name = "offset")
quadratic_coeff <- Variable(1, name = "quadratic_coeff")

Tlen = 100
position <- Variable(2, Tlen, name = "position")
velocity <- Variable(2, Tlen, name = "velocity")
force <- Variable(2, Tlen - 1, name = "force")

xs <- Variable(80, name = "xs")
xsr <- Variable(200, name = "xsr")
xef <- Variable(80, name = "xef")

# Check for all installed QP solvers
solvers <- installed_solvers()
solvers <- solvers[solvers %in% qp_solvers()]
if("MOSEK" %in% installed_solvers())
  solvers <- c(solvers, "MOSEK")

solve_QP <- function(problem, solver_name) {
  solve(problem, solver = solver_name, verbose = TRUE)
}

test_quad_over_lin <- function(solver) {
  p <- Problem(Minimize(0.5 * quad_over_lin(abs(x-1), 1)), list(x <= -1))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), as.matrix(c(-1, -1)), tolerance = 1e-4)
  
  for(con in constraints(p))
    expect_equal(result$getDualValue(con), as.matrix(c(2, 2)), tolerance = 1e-4)
}

test_abs <- function(solver) {
  u <- Variable(2)
  constr <- list()
  constr <- c(constr, abs(u[2] - u[1]) <= 100)
  prob <- Problem(Minimize(sum_squares(u)), constr)
  print(paste("The problem is QP: ", is_qp(prob), sep = ""))
  expect_true(is_qp(prob))
  result <- solve(prob, solver = solver)
  expect_equal(result$value, 0, tolerance = TOL)
}

test_power <- function(solver) {
  p <- Problem(Minimize(sum(power(x, 2))), list())
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), as.matrix(c(0, 0)), tolerance = 1e-4)
}

test_power_matrix <- function(solver) {
  p <- Problem(Minimize(sum(power(A - 3., 2))), list())
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), as.matrix(rep(3, 4)), tolerance = 1e-4)
}

test_square_affine <- function(solver) {
  A <- matrix(rnorm(10*2), nrow = 10, ncol = 2)
  b <- matrix(rnorm(10), nrow = 10, ncol = 1)
  p <- Problem(Minimize(sum_squares(A %*% x - b)))
  result <- solve_QP(p, solver)
  
  x_star <- base::solve(qr(A),b)
  for(var in variables(p))
    expect_equal(result$getValue(var), x_star, tolerance = 0.1)
}

test_quad_form <- function(solver) {
  A <- matrix(rnorm(5*5), nrow = 5, ncol = 5)
  z <- matrix(rnorm(5), nrow = 5, ncol = 1)
  P <- t(A) %*% A
  q <- -2*P %*% z
  p <- Problem(Minimize(quad_form(w, P) + t(q) %*% w))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), z, tolerance = 1e-4)
}

test_affine_problem <- function(solver) {
  A <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  A <- pmax(A, 0)
  b <- matrix(rnorm(5), nrow = 5, ncol = 1)
  b <- pmax(b, 0)
  p <- Problem(Minimize(sum(x)), list(x >= 0, A %*% x <= b))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), as.matrix(c(0,0)), tolerance = 1e-3)
}

test_maximize_problem <- function(solver) {
  A <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  A <- pmax(A, 0)
  b <- matrix(rnorm(5), nrow = 5, ncol = 1)
  b <- pmax(b, 0)
  p <- Problem(Maximize(-sum(x)), list(x >= 0, A %*% x <= b))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), as.matrix(c(0,0)), tolerance = 1e-3)
}

test_norm_2 <- function(solver) {
  A <- matrix(rnorm(10*5), nrow = 10, ncol = 5)
  b <- matrix(rnorm(10), nrow = 10)
  p <- Problem(Minimize(p_norm(A %*% w - b, 2)))
  result <- solve_QP(p, solver)
  
  x_star <- base::solve(qr(A),b)
  for(var in variables(p))
    expect_equal(result$getValue(var), x_star, tolerance = 0.1)
}

test_mat_norm_2 <- function(solver) {
  A <- matrix(rnorm(5*3), nrow = 5, ncol = 3)
  B <- matrix(rnorm(5*2), nrow = 5, ncol = 2)
  p <- Problem(Minimize(p_norm(A %*% C - B, 2)))
  result <- solve_QP(p, solver)
  
  C_star <- base::solve(qr(A),B)
  for(var in variables(p))
    expect_equal(result$getValue(var), C_star, tolerance = 0.1)
}

test_quad_form_coeff <- function(solver) {
  A <- matrix(rnorm(5*5), nrow = 5, ncol = 5)
  z <- matrix(rnorm(5), nrow = 5, ncol = 1)
  P <- t(A) %*% A
  q <- -2*P %*% z
  p <- Problem(Minimize(quad_form(w, P) + t(q) %*% w))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), z, tolerance = 1e-4)
}

test_quad_form_bound <- function(solver) {
  P <- rbind(c(13, 12, -2), c(12, 17, 6), c(-2, 6, 12))
  q <- as.matrix(c(-22, -14.5, 13))
  r <- 1
  y_star <- as.matrix(c(1, 0.5, -1))
  p <- Problem(Minimize(0.5*quad_form(y, P) + t(q) %*% y + r),
               list(y >= -1, y <= 1))
  result <- solve_QP(p, solver)
  for(var in variables(p))
    expect_equal(result$getValue(var), y_star, tolerance = 1e-4)
}

test_that("test all solvers", {
  for(solver in solvers) {
    test_quad_over_lin(solver)
    test_power(solver)
    test_power_matrix(solver)
    test_square_affine(solver)
    test_quad_form(solver)
    test_affine_problem(solver)
    test_maximize_problem(solver)
    test_abs(solver)
    
    # Do we need the following functionality?
    # test_norm_2(solver)
    # test_mat_norm_2(solver)
    
    test_quad_form_coeff(solver)
    test_quad_form_bound(solver)
    test_regression_1(solver)
    test_regression_2(solver)
    
    # Slow tests:
    test_control(solver)
    test_sparse_system(solver)
    test_smooth_ridge(solver)
    test_huber_small(solver)
    test_huber(solver)
    test_equivalent_forms_1(solver)
    test_equivalent_forms_2(solver)
    test_equivalent_forms_3(solver)
  }
})
