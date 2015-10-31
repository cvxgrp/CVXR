test_that("Find the largest Euclidean ball in the polyhedron", {
  # Create the input data
  a1 <- matrix(c(2,1))
  a2 <- matrix(c(2,-1))
  a3 <- matrix(c(-1,2))
  a4 <- matrix(c(-1,-2))
  b <- rep(1,4)
  
  # Create and solve the model
  r <- Variable(name = "r")
  x_c <- Variable(2, name = "x_c")
  obj <- Maximize(r)
  constraints <- list(
    t(a1)*x_c + norm(a1,"F")*r <= b[1],
    t(a2)*x_c + norm(a2,"F")*r <= b[2],
    t(a3)*x_c + norm(a3,"F")*r <= b[3],
    t(a4)*x_c + norm(a4,"F")*r <= b[4]
  )
  
  p <- Problem(obj, constraints)
})

test_that("Test examples from the README", {
  # Problem data
  m <- 30
  n <- 20
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m, ncol = 1)
  
  # Construct the problem
  x <- Variable(n)
  objective <- Minimize(SumEntries(Square(A*x - b)))
  constraints <- list(x >= 0, x <= 1)
  p <- Problem(objective, constraints)

  ###########################################
  # Scalar variable
  a <- Variable()
  
  # Column vector variable of length 5
  x <- Variable(5)
  
  # Matrix variable with 4 rows and 7 columns
  A <- Variable(4, 7)
  ###########################################
  # Positive scalar parameter
  m <- Parameter(sign = "positive")
  
  # Column vector parameter with unknown sign (by default)
  c <- Parameter(5)
  
  # Matrix parameter with negative entries
  G <- Parameter(4, 7, sign = "negative")
  
  ###########################################
  a <- Variable()
  x <- Variable(5)
  
  # expr is an Expression object after each assignment
  expr <- 2*x
  expr <- expr - a
  expr <- SumEntries(expr) + norm(x, 2)
  ###########################################
  # Problem data
  n <- 10
  m <- 5
  A <- matrix(rnorm(n*m), nrow = n, ncol = m)
  b <- matrix(rnorm(n), nrow = n, ncol = 1)
  gamma <- Parameter(sign = "positive")
  
  # Construct the problem
  x <- Variable(m)
  objective <- Minimize(SumEntries(Square(A*x - b)) + gamma*norm(x, 1))
  p <- Problem(objective)
  ###########################################
  n <- 10
  
  mu <- matrix(rnorm(n), nrow = 1, ncol = n)
  sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
  sigma <- t(sigma) * sigma
  gamma <- Parameter(sign = "positive")
  gamma@value <- 1
  x <- Variable(n)
  
  expected_return <- mu*x
  # risk <- QuadForm(x, sigma)
  
  # objective <- Maximize(expected_return - gamma*risk)
  # p <- Problem(objective, list(SumEntries(x) == 1))
  ###########################################
  N <- 50
  M <- 40
  n <- 10
  data1 <- lapply(1:N, function(i) { list(1, matrix(rnorm(n, mean = 1.0, sd = 2.0))) })
  data2 <- lapply(1:M, function(i) { list(-1, matrix(rnorm(n, mean = -1.0, sd = 2.0))) })
  data <- c(data1, data2)
  
  # Construct problem
  gamma <- Parameter(sign = "positive")
  gamma@value <- 0.1
  a <- Variable(n)
  b <- Variable()
  
  slack <- lapply(data, function(x) { 
    label <- x[[1]]
    sample <- x[[2]]
    Pos(1 - label*(t(sample)*a - b))
  })
  # objective <- Minimize(norm(a, 2) * gamma * Reduce("+", slack))
  # p <- Problem(objective)
})


