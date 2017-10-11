TOL <- 1e-6

test_that("Test quadratic form with a singular matrix", {
  require(Matrix)
  
  # Solve a quadratic program
  for(n in c(3,4,5)) {
    for(i in 0:4) {
      # Construct a random 1-d finite distribution
      v <- exp(rnorm(n))
      v <- v / sum(v)
      
      # Construct a random positive definite matrix
      A <- matrix(rnorm(n^2), nrow = n, ncol = n)
      Q <- A %*% t(A)
      
      # Project onto the orthogonal complement of v
      # This turns Q into a singular matrix with a known nullspace
      E <- diag(rep(1,n)) - v %*% t(v) / as.numeric(t(v) %*% v)
      Q <- E %*% (Q %*% t(E))
      Q_rank <- rankMatrix(Q)
      observed_rank <- Q_rank[1]
      desired_rank <- n-1
      expect_equal(observed_rank, desired_rank)
      
      for(action in c("minimize", "maximize")) {
        # Look for the extremum of the quadratic form under the simplex constraint
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
        
        # Check that CVXR found the right answer
        xopt <- result$getValue(x)
        yopt <- t(xopt) %*% (Q %*% xopt)
        sapply(yopt, function(y) { expect_equal(y, 0, tolerance = 1e-3) })
        mapply(function(x_elem, v_elem) { expect_equal(x_elem, v_elem, tolerance = 1e-3) }, xopt, v)
      }
    }
  }
})

test_that("Test quadratic form with a sparse matrix", {
  require(Matrix)
  Q <- sparseMatrix(i = 1:2, j = 1:2, x = rep(1, 2))
  x <- Variable(2)
  cost <- quad_form(x, Q)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob)
  expect_equal(result$value, 5, tolerance = TOL)
})

test_that("Test when P is constant and not symmetric", {
  P <- rbind(c(2, 2), c(3, 4))
  x <- Variable(2)
  cost <- quad_form(x, P)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  result <- solve(prob)
  expect_equal(result$value, 28, tolerance = TOL)
})

test_that("Test error when P is symmetric but not definite", {
  P <- rbind(c(1, 0), c(0, -1))
  x <- Variable(2)
  
  # Forming quadratic form is okay
  expect_warning(cost <- quad_form(x, P))
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  expect_error(solve(prob))
})
