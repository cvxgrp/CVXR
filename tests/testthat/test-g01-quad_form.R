context("test-g01-quad_form")
TOL <- 1e-6

test_that("Test quadratic form with a singular matrix", {
  skip_on_cran()
  ## Solve a quadratic program
  for(n in c(3,4,5)) {
      for(i in 0:4) {
          ## Construct a random 1-d finite distribution
          v <- exp(stats::rnorm(n))
          v <- v / sum(v)

          ## Construct a random positive definite matrix
          A <- matrix(stats::rnorm(n^2), nrow = n, ncol = n)
          Q <- A %*% t(A)

          ## Project onto the orthogonal complement of v
          ## This turns Q into a singular matrix with a known nullspace
          E <- diag(rep(1,n)) - v %*% t(v) / as.numeric(t(v) %*% v)
          Q <- E %*% (Q %*% t(E))
          Q_rank <- Matrix::rankMatrix(Q)
          observed_rank <- Q_rank[1]
          desired_rank <- n-1
          expect_equal(observed_rank, desired_rank)

          for(action in c("minimize", "maximize")) {
              ## Look for the extremum of the quadratic form under the simplex constraint
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
              result <- solve(p)

              ## Check that CVXR found the right answer
              xopt <- result$getValue(x)
              yopt <- t(xopt) %*% (Q %*% xopt)
              sapply(yopt, function(y) { expect_equal(y, 0, tolerance = 1e-3) })
              mapply(function(x_elem, v_elem) { expect_equal(x_elem, v_elem, tolerance = 1e-3) }, xopt, v)
          }
      }
  }
})

test_that("Test quadratic form with a sparse matrix", {
  skip_on_cran()
  Q <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = rep(1, 2))
  x <- Variable(2)
  cost <- quad_form(x, Q)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob)
  expect_equal(result$value, 5, tolerance = TOL)

  # Here are our QP factors
  A <- Constant(Matrix::sparseMatrix(i = 1:4, j = 1:4, x = rep(1, 4)))
  c <- matrix(1, nrow = 1, ncol = 4)

  # Here is our optimization variable
  x <- Variable(4)

  # And the QP problem setup
  fun <- quad_form(x, A) - c %*% x
  objective <- Minimize(fun)
  problem <- Problem(objective)

  result <- solve(problem)
  expect_equal(length(result$getValue(fun)), 1)
})

# test_that("Test quadratic form with a parameter", {
#   P <- Parameter(2,2,PSD = TRUE)
#   Q <- diag(2)
#   x <- Variable(2)
#   cost <- quad_form(x, P)
#   value(P) <- Q
#   prob <- Problem(Minimize(cost), list(x == c(1,2)))
#   # result <- solve(prob)
#   # expect_equal(result$value, 5, tolerance = TOL)
# })

test_that("Test when P is constant and not symmetric", {
  skip_on_cran()
  P <- rbind(c(2, 2), c(3, 4))
  x <- Variable(2)
  cost <- quad_form(x, P)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  expect_error(solve(prob), "Problem does not follow DCP rules.")
})

test_that("Test error when P is symmetric but not definite", {
  skip_on_cran()
  P <- rbind(c(1, 0), c(0, -1))
  x <- Variable(2)

  # Forming quadratic form is okay
  cost <- quad_form(x, P)
  # expect_warning(cost <- quad_form(x, P))
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  expect_error(solve(prob), "Problem does not follow DCP rules.")
})

test_that("Test case where objective evaluation differs from result", {
  skip_on_cran()
  x <- Variable(2,1)
  A <- matrix(1, nrow = 1, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 1)
  obj0 <- -t(B) %*% x
  obj1 <- quad_form(t(B) %*% x, A)
  prob <- Problem(Minimize(obj0 + obj1))
  result <- solve(prob)
  expect_equal(result$value, result$getValue(prob@objective@expr))
})

test_that("Test a quadratic form multiplied by zero", {
  skip_on_cran()
  data_norm <- runif(5)
  M <- matrix(runif(5*2), nrow = 5, ncol = 2)
  c <- Variable(ncol(M))
  lopt <- 0
  laplacian_matrix <- matrix(1, nrow = 2, ncol = 2)
  design_matrix <- Constant(M)
  objective <- Minimize(sum_squares(design_matrix %*% c - data_norm) + lopt*quad_form(c, laplacian_matrix))
  constraints <- list((M[1,1] * c) == 1)   # (K*c >= -0.1)
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  expect_equal(result$status, "optimal")
})
