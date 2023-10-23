context("test-g01-quad_form")
TOL <- 1e-6
EIGVAL_TOL <- CVXR:::EIGVAL_TOL

test_that("Test quad_form with a singular matrix", {
  skip_on_cran()
  
  # Solve a quadratic program.
  set.seed(1234)
  for(n in c(3,4,5)) {
      for(i in 0:4) {
          # Construct a random 1-D finite distribution.
          v <- exp(stats::rnorm(n))
          v <- v / sum(v)

          # Construct a random positive definite matrix.
          A <- matrix(stats::rnorm(n^2), nrow = n, ncol = n)
          Q <- A %*% t(A)

          # Project onto the orthogonal complement of v.
          # This turns Q into a singular matrix with a known nullspace.
          E <- diag(n) - v %*% t(v) / as.numeric(t(v) %*% v)
          Q <- E %*% (Q %*% t(E))
          Q_rank <- Matrix::rankMatrix(Q)
          observed_rank <- Q_rank[1]
          desired_rank <- n - 1
          expect_equal(observed_rank, desired_rank)

          for(action in c("minimize", "maximize")) {
              # Look for the extremum of the quadratic form under the simplex constraint.
              x <- Variable(n)
              if(action == "minimize") {
                  q <- quad_form(x, Q)
                  objective <- Minimize(q)
              } else if(action == "maximize") {
                  q <- quad_form(x, -Q)
                  objective <- Maximize(q)
              }
              constraints <- list(x >= 0, sum(x) == 1)
              p <- Problem(objective, constraints)
              result <- solve(p, solver = "OSQP")

              # Check that CVXR found the right answer.
              xopt <- result$getValue(x)
              yopt <- t(xopt) %*% (Q %*% xopt)
              expect_true(is.allclose(yopt, 0, atol = 1e-3))
              expect_true(is.allclose(xopt, v, atol =))
          }
      }
  }
})

test_that("Test quad_form with a sparse matrix", {
  skip_on_cran()
  Q <- Matrix(diag(2), sparse = TRUE)
  x <- Variable(2)
  cost <- quad_form(x, Q)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob, solver = "OSQP")
  expect_equal(result$value, 5, tolerance = TOL)

  # Here are our QP factors.
  A <- Constant(Matrix(diag(4), sparse = TRUE))
  c <- matrix(1, nrow = 1, ncol = 4)

  # Here is our optimization variable.
  x <- Variable(4)

  # And the QP problem setup.
  fun <- quad_form(x, A) - c %*% x
  objective <- Minimize(fun)
  problem <- Problem(objective)

  result <- solve(problem, solver = "OSQP")
  expect_equal(length(result$getValue(fun)), 1)
})

test_that("Test quad_form with a parameter", {
  skip_on_cran()
  P <- Parameter(2,2,PSD = TRUE)
  Q <- diag(2)
  x <- Variable(2)
  cost <- quad_form(x, P)
  value(P) <- Q
  prob <- Problem(Minimize(cost), list(x == c(1,2)))
  result <- solve(prob, solver = "SCS")
  # expect_warning(result <- solve(prob, solver = "SCS"))
  expect_equal(result$value, 5, tolerance = TOL)
})

test_that("Test when P is constant and not symmetric", {
  skip_on_cran()
  P <- rbind(c(2, 2), c(3, 4))
  x <- Variable(2)
  expect_error(quad_form(x, P), "Quadratic form matrices must be symmetric/Hermitian.", fixed = TRUE)
})

test_that("Test error when P is symmetric but not definite", {
  skip_on_cran()
  P <- rbind(c(1, 0), c(0, -1))
  x <- Variable(2)

  # Forming quad_form is okay
  cost <- quad_form(x, P)
  # expect_warning(cost <- quad_form(x, P))
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  expect_error(solve(prob, solver = "SCS"), "Problem does not follow DCP rules.", fixed = TRUE)
})

test_that("Test that PSD check when eigenvalue is exactly -EIGVAL_TOL", {
  skip_on_cran()
  P <- rbind(c(-0.999*EIGVAL_TOL, 0), c(0, 10))
  x <- Variable(2)
  
  # Forming quad_form is okay.
  cost <- quad_form(x, P)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob, solver = "SCS")
})

test_that("Test that NSD check when eigenvalue is exactly EIGVAL_TOL", {
  skip_on_cran()
  P <- rbind(c(0.999*EIGVAL_TOL, 0), c(0, -10))
  x <- Variable(2)
  
  # Forming quad_form is okay.
  cost <- quad_form(x, P)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob, solver = "SCS")
})

test_that("Test case where objective evaluation differs from result", {
  skip_on_cran()
  x <- Variable(2,1)
  A <- matrix(1.0)
  B <- matrix(1.0, nrow = 2, ncol = 1)
  obj0 <- -t(B) %*% x
  obj1 <- quad_form(t(B) %*% x, A)
  prob <- Problem(Minimize(obj0 + obj1))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, result$getValue(prob@objective@expr))
})

test_that("Test a quad_form multiplied by zero", {
  skip_on_cran()
  data_norm <- runif(5)
  M <- matrix(runif(5*2), nrow = 5, ncol = 2)
  c <- Variable(ncol(M))
  lopt <- 0
  laplacian_matrix <- matrix(1, nrow = 2, ncol = 2)
  design_matrix <- Constant(M)
  objective <- Minimize(sum_squares(design_matrix %*% c - data_norm) + lopt*quad_form(c, laplacian_matrix))
  constraints <- list((M[1] %*% c) == 1)   # (K*c) >= -0.1
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "SCS")
})

test_that("Test quad_form with P = 0", {
  skip_on_cran()
  x <- Variable(3)
  A <- diag(3)
  b <- matrix(rep(1,3))
  c <- -matrix(rep(1,3))
  P <- matrix(0, nrow = 3, ncol = 3)
  expr <- (1/2)*quad_form(x, P) + t(c) %*% x
  prob <- Problem(Minimize(expr), list(A %*% x <= b))
  solve(prob, solver = "SCS")
})

test_that("test assume_PSD argument", {
  skip_on_cran()
  x <- Variable(3)
  A <- diag(3)
  expr <- quad_form(x, A, assume_PSD = TRUE)
  expect_true(is_convex(expr))
  
  A <- -diag(3)
  expr <- quad_form(x, A, assume_PSD = TRUE)
  expect_true(is_convex(expr))
  
  prob <- Problem(Minimize(expr))
  # Transform to a SolverError.
  expect_error(solve(prob, solver = "OSQP"), "Workspace allocation error!")
})
