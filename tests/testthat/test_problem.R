a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the variables method", {
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars_ <- variables(p)
  ref <- list(a, x, b, A)
  mapply(function(v, r) { expect_equal(v, r) }, vars_, ref)
})

test_that("test the parameters method", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(p)
  ref <- c(p1, p2, p3)
  mapply(function(p, r) { expect_equal(p, r) }, params, ref)
})

test_that("test the constants method", {
  c1 <- matrix(randn(2), nrow = 1, ncol = 2)
  c2 <- matrix(randn(2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(c1*x), list(x >= c2))
  # constants_ <- constants(p)
  # ref <- list(as.character(c1), as.character(c2))
  # mapply(function(c, r) { expect_equal(c, r) }, constants_, ref)
})

test_that("Test the size_metrics method", {
  
})

test_that("Test the get_problem_data method", {
  
})

test_that("Test the unpack_results method", {
  
})

test_that("Test silencing and enabling solver messages", {
  
})

test_that("Test registering other solve methods", {
  
})

test_that("Test that variables and constraints keep a consistent order", {
  
})

test_that("Test removing duplicate constraints objects", {
  
})

test_that("test the is_dcp method", {
  p <- Problem(Minimize(NormInf(a)))
  expect_true(is_dcp(p))
  
  p <- Problem(Maximize(NormInf(a)))
  expect_false(is_dcp(p))
  # expect_error(solve(p))
  # solve(p, ignore_dcp = TRUE)
})

test_that("test problems involving variables with the same name", {
  var <- Variable(name = "a")
  p <- Problem(Maximize(a + var), list(var == 2 + a, var <= 3))
  # result <- solve(p)
  # expect_equal(result, 4.0, tolerance = TOL)
  # expect_equal(result$a, 1, tolerance = TOL)
  # expect_equal(result$var, 3, tolerance = TOL)
})

test_that("test adding problems", {
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob_minimize <- prob1 + prob2
  expect_equal(length(prob_minimize@constraints), 3)
  # result <- solve(prob_minimize)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)

  prob3 <- Problem(Maximize(a), list(b <= 1))
  prob4 <- Problem(Maximize(2*b), list(a <= 2))
  prob_maximize <- prob3 + prob4
  expect_equal(length(prob_maximize@constraints), 2)
  # result <- solve(prob_maximize)
  # expect_equal(result$optimal_value, 4, tolerance = TOL)
  
  # Test using the sum function
  prob5 <- Problem(Minimize(3*a))
  prob_sum <- Reduce("+", list(prob1, prob2, prob5))
  expect_equal(length(prob_sum@constraints), 3)
  # result <- solve(prob_sum)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  prob_sum <- Reduce("+", list(prob1))
  expect_equal(length(prob_sum@constraints), 1)
  
  # Test Minimize + Maximize
  expect_error(prob_bad_sum <- prob1 + prob3)
})

test_that("test problem multiplication by scalar", {
  prob1 <- Problem(Minimize(a^2), list(a >= 2))
  # answer <- solve(prob1)
  factors <- c(0, 1, 2.3, -4.321)
  for(f in factors) {
    # expect_equal(solve(f * prob1)$optimal_value, f * answer, tolerance = TOL)
    # expect_equal(solve(prob1 * f)$optimal_value, f * answer, tolerance = TOL)
  }
})

test_that("test problem linear combinations", {
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob3 <- Problem(Maximize(-(b + a)^2), list(b >= 3))
  
  # Simple addition and multiplication
  combo1 <- prob1 + 2 * prob2
  combo1_ref <- Problem(Minimize(a + 4*b), list(a >= b, a >= 1, b >= 2))
  # expect_equal(solve(combo1), solve(combo1_ref))
  
  # Division and subtraction
  combo2 <- prob1 - prob3/2
  combo2_ref <- Problem(Minimize(a + (b + a)^2), list(b >= 3, a >= b))
  # expect_equal(solve(combo2), solve(combo2_ref))
  
  # Multiplication with 0 (prob2's constraints should still hold)
  combo3 <- prob1 + 0*prob2 - 3*prob3
  combo3_ref <- Problem(Minimize(a + 3*(b + a)^2), list(a >= b, a >= 1, b >= 3))
  # expect_equal(solve(combo3), solve(combo3_ref))
})

test_that("test solving problems in parallel", {
  p <- Parameter()
  problem <- Problem(Minimize(Square(a) + Square(b) + p), list(b >= 2, a >= 1))
  p@value <- 1
  
  # Ensure that parallel solver still works after repeated calls
  for(i in 1:2) {
    # result <- solve(problem, parallel = TRUE)
    # expect_equal(result$optimal_value, 6.0, tolerance = TOL)
    # expect_equal(result$status, "OPTIMAL")
    # expect_equal(result$a, 1, tolerance = TOL)
    # expect_equal(result$b, 2, tolerance = TOL)
  }
  
  # The constant p should not be a separate problem, but rather added to the first separable problem
  # expect_true(length(problem@.separable_problems) == 2)
  
  # Ensure that parallel solver works with options
  # result <- solve(problem, parallel = TRUE, verbose = TRUE, warm_start = TRUE)
  # expect_equal(result$optimal_value, 6.0, tolerance = TOL)
  # expect_equal(result$status, "OPTIMAL")
  # expect_equal(result$a, 1, tolerance = TOL)
  # expect_equal(result$b, 2, tolerance = TOL)
  
  # Ensure that parallel solver works when problem changes
  problem@objective <- Minimize(Square(a) + Square(b))
  # result <- solve(problem, parallel = TRUE)
  # expect_equal(result$optimal_value, 5.0, tolerance = TOL)
  # expect_equal(result$status, "OPTIMAL")
  # expect_equal(result$a, 1, tolerance = TOL)
  # expect_equal(result$b, 2, tolerance = TOL)
})

test_that("Test scalar LP problems", {
  p <- Problem(Minimize(3*a), list(a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  
  p <- Problem(Maximize(3*a - b), list(a <= 2, b == a, b <= 5))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 4.0, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  # expect_equal(result$b, 2, tolerance = TOL)
  
  # With a constant in the objective
  p <- Problem(Minimize(3*a - b + 100), list(a >= 2, b + 5*c - 2 == a, b <= 5 + c))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 101+1.0/6, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  # expect_equal(result$b, 5-1.0/6, tolerance = TOL)
  # expect_equal(result$c, -1.0/6, tolerance = TOL)
  
  # Test status and value
  exp <- Maximize(a)
  p <- Problem(exp, list(a <= 2))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, p@value)
  # expect_equal(result$status, "OPTIMAL")
  # expect_false(is.na(result$a))
  # expect_false(is.na(p@constraints[1]@dual_value))
  
  # Unbounded problems
  p <- Problem(Maximize(a), list(a >= 2))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, Inf)
  # expect_true(result$optimal_value > 0)
  # expect_true(is.na(result$a))
  # expect_true(is.na(p@constraints[1]@dual_value))
  
  # Infeasible problems
  p <- Problem(Maximize(a), list(a >= 2, a <= 1))
  # a.save_value(2)
  # p.constraints[1].save_value(2)
  
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, p@value)
  # expect_equal(result$status, "INFEASIBLE")
  # expect_equal(result$optimal_value, Inf)
  # expect_true(result$optimal_value < 0)
  # expect_true(is.na(result$a))
  # expect_true(is.na(p@constraints[1]@dual_value))
  
  p <- Problem(Minimize(-a), list(a >= 2, a <= 1))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, p@value)
  # expect_equal(result$status, "INFEASIBLE")
  # expect_true(result$optimal_value, Inf)
  # expect_true(result$optimal_value > 0)
})

test_that("Test vector LP problems", {
  c <- Constant(matrix(c(1, 2), nrow = 2, ncol = 1))@value
  p <- Problem(Minimize(t(c) * x), list(x >= c))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 5, tolerance = TOL)
  # expect_equal(result$x, c(1, 2), tolerance = TOL)
  
  A <- Constant(rbind(c(3, 5), c(1, 2)))@value
  I <- Constant(rbind(c(1, 0), c(0, 1)))
  p <- Problem(Minimize(t(c) * x + a), list(A*x >= c(-1, 1), 4*I*z == x, z >= c(2,2), a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 26, tolerance = 1e-3)
  # obj <- t(c) * result$x + result$a
  # expect_equal(obj[1,1], result$optimal_value, tolerance = TOL)
  # expect_equal(result$x, c(8,8), tolerance = 1e-3)
  # expect_equal(result$z, c(2,2), tolerance = 1e-3)
})

test_that("Test ECOS with no inequality constraints", {
  T <- Constant(matrix(1, nrow = 2, ncol = 2))@value
  p <- Problem(Minimize(1), list(A == T))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  # expect_equal(result$A, T, tolerance = TOL)
})

test_that("Test matrix LP problems", {
  T <- Constant(matrix(1, nrow = 2, ncol = 2))@value
  p <- Problem(Minimize(1), list(A == T))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  # expect_equal(result$A, T, tolerance = TOL)
  
  T <- Constant(matrix(1, nrow = 2, ncol = 3)*2)@value
  c <- Constant(matrix(c(3,4), nrow = 2, ncol = 1))@value
  p <- Problem(Minimize(1), list(A >= T*C, A == B, C == t(T)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  # expect_equal(result$A, result$B, tolerance = TOL)
  # expect_equal(result$C, T, tolerance = TOL)
  # expect_true(all(result$A >= T*result$C))
  
  # Test variables are dense
  # expect_true(is(A, "matrix"))
  
})

test_that("Test variable promotion", {
  p <- Problem(Minimize(a), list(x <= a, x == c(1, 2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  
  p <- Problem(Minimize(a), list(A <= a, A == rbind(c(1,2), c(3,4))))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 4, tolerance = TOL)
  # expect_equal(result$a, 4, tolerance = TOL)
  
  # Promotion must happen before multiplication
  p <- Problem(Minimize(matrix(1, nrow = 1, ncol = 2)*(x + a + 1)), list(a + x >= c(1, 2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 5, tolerance = TOL)
})

test_that("Test parameter promotion", {
  a <- Parameter()
  exp <- rbind(c(1,2), c(3,4))*a
  a@value <- 2
  expect_false(any(exp@value - 2*cbind(c(1,2), c(3,4))))
})

test_that("test problems with parameters", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  p1@value <- 2
  p2@value <- -matrix(1, nrow = 3, ncol = 1)
  p3@value <- matrix(1, nrow = 4, ncol = 4)
  result <- solve(p)
  expect_equal(result$optimal_value, -6, tolerance = TOL)
  
  p1@value <- NA
  expect_error(solve(p))
})

test_that("test problems with NormInf", {
  # Constant argument
  p <- Problem(Minimize(NormInf(-2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  # Scalar arguments
  p <- Problem(Minimize(NormInf(a)), list(a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
  
  p <- Problem(Minimize(3*NormInf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 3, tolerance = TOL)
  # expect_equal(result$a + 2*result$b, 0, tolerance = TOL)
  # expect_equal(result$c, 3, tolerance = TOL)

  # Maximize
  p <- Problem(Maximize(-NormInf(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, -2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Vector arguments
  p <- Problem(Minimize(NormInf(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  # expect_equal(result$x[2] - result$z[2], 7, tolerance = TOL)
})

test_that("Test problems with Norm1", {
  # Constant argument
  p <- Problem(Minimize(Norm1(-2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  # Scalar arguments
  p <- Problem(Minimize(Norm1(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Maximize
  p <- Problem(Maximize(-Norm1(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, -2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Vector arguments
  p <- Problem(Minimize(Norm1(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 15, tolerance = TOL)
  # expect_equal(result$x[2] - result$z[2], 7, tolerance = TOL)
})

test_that("Test problems with Norm2", {
  # Constant argument
  p <- Problem(Minimize(Norm2(-2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  # Scalar arguments
  p <- Problem(Minimize(Norm2(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Maximize
  p <- Problem(Maximize(-Norm2(a)), list(a <= -2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, -2, tolerance = TOL)
  # expect_equal(result$a, -2, tolerance = TOL)
  
  # Vector arguments
  p <- Problem(Minimize(Norm2(x - z) + 5), list(x <= c(2,3), z <= c(-1,-4)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 12.61577, tolerance = TOL)
  # expect_equal(result$x, c(2,3), tolerance = TOL)
  # expect_equal(result$z, c(-1,-4), tolerance = TOL)
  
  # Row arguments
  p <- Problem(Minimize(Norm2(t(x - z)) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 12.61577, tolerance = TOL)
  # expect_equal(result$x, c(2,3))
  # expect_equal(result$z, c(-1,-4))
})

test_that("Test problems with abs", {
  p <- Problem(Minimize(SumEntries(abs(A))), list(-2 >= A))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 8, tolerance = TOL)
  # expect_equal(result$A, rep(-2, 4), tolerance = TOL)
})

test_that("Test problems with QuadForm", {
  # expect_error(solve(Problem(Minimize(QuadForm(x, A)))))
  # expect_error(solve(Problem(Minimize(QuadForm(1, A)))))
  # expect_error(solve(Problem(Minimize(QuadForm(x, rbind(c(-1,0), c(0.9)))))))
  
  P <- rbind(c(4,0), c(0,9))
  p <- Problem(Minimize(QuadForm(x, P)), list(x >= 1))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 13, tolerance = 1e-3)
  
  c <- c(1,2)
  p <- Problem(Minimize(QuadForm(c, A)), list(A >= 1))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 9, tolerance = TOL)
  
  c <- c(1,2)
  P <- rbind(c(4,0), c(0,9))
  p <- Problem(Minimize(QuadForm(c, P)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 40)
})

test_that("Test combining atoms", {
  p <- Problem(Minimize(Norm2(5 + Norm1(z) + Norm1(x) + NormInf(x - z))), list(x >= c(2,3), z <= c(-1,-4), Norm2(x + z) <= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 22, tolerance = TOL)
  # expect_equal(result$x, c(2,3))
  # expect_equal(result$z, c(-1,-4))
})

test_that("Test multiplying by constant atoms", {
  p <- Problem(Minimize(Norm2(c(3,4)) * a), list(a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(result$a, 2, tolerance = TOL)
})

test_that("Test recovery of dual variables", {
  for(solver in c("ECOS", "SCS")) {
    if(solver == "SCS")
      acc <- 1e-1
    else
      acc <- 1e-5
    P <- Problem(Minimize(Norm1(x + z)), list(x >= c(2,3), rbind(c(1,2), c(3,4)) * z == c(-1,4), Norm2(x + z) <= 100))
    result <- solve(p, solver = solver)
    expect_equal(result$optimal_value, 4, tolerance = acc)
    expect_equal(result$x, c(4,3), tolerance = acc)
    expect_equal(result$z, c(-4,1), tolerance = acc)
    
    # Dual values
    expect_equal(result$constraints[1]$dual_value, c(0,1), tolerance = acc)
    expect_equal(result$constraints[2]$dual_value, c(-1,0.5), tolerance = acc)
    expect_equal(result$constraints[3]$dual_value, 0, tolerance = acc)
    
    Tmat <- matrix(1, nrow = 2, ncol = 3) * 2
    c <- matrix(c(3,4), nrow = 1, ncol = 2)
    p <- Problem(Minimize(1), list(A >= Tmat*C, A == B, C == t(Tmat)))
    result <- solve(p, solver = solver)
    
    # Dual values
    expect_equal(result$constraints[1]$dual_value, rep(0,4), tolerance = acc)
    expect_equal(result$constraints[2]$dual_value, rep(0,4), tolerance = acc)
    expect_equal(result$constraints[3]$dual_value, rep(0,6), tolerance = acc)
  }
})

test_that("Test problems with indexing", {
  # Vector variables
  p <- Problem(Maximize(x[1,1]), list(x[1,1] <= 2, x[2,1] == 3))
  result <- solve(p)
  expect_equal(result$optimal_value, 2, tolerance = TOL)
  expect_equal(result$x, c(2,3))
  
  n <- 10
  A <- matrix(0:(n^2-1), nrow = n, ncol = n)
  x <- Variable(n,n)
  p <- Problem(Minimize(SumEntries(x)), list(x == A))
  result <- solve(p)
  answer <- n*n*(n*n+1)/2 - n*n
  expect_equal(result$optimal_value, answer)
  
  # Matrix variables
  obj <- A[1,1] + A[1,2] + A[2,2] + A[2,1]
  p <- Problem(Maximize(obj), list(A <= rbind(c(1,-2), c(-3,4))))
  result <- solve(p)
  expect_equal(result$optimal_value, 0, tolerance = TOL)
  expect_equal(result$A, c(1,-2,-3,4))
  
  # Indexing arithmetic expressions
  exp <- rbind(c(1,2), c(3,4)) * z + x
  p <- Problem(Minimize(exp[2,1]), list(x == z, z == c(1,2)))
  result <- solve(p)
  expect_equal(result$optimal_value, 12, tolerance = TOL)
  expect_equal(result$x, result$z, tolerance = TOL)
})

test_that("Test problems with slicing", {
  p <- Problem(Maximize(SumEntries(C)), list(C[2:4,] <= 2, C[1,] == 1))
  result <- solve(p)
  expect_equal(result$optimal_value, 10, tolerance = TOL)
  expect_equal(result$C, 2*c(1,2,2))
  
  # TODO: Finish this
})

test_that("Test using a CVXR expression as an objective", {
  expect_error(Problem(x+2))
})

test_that("Test multiplication by zero", {
  exp <- 0*a
  expect_equal(value(exp), 0)
  obj <- Minimize(exp)
  p <- Problem(obj)
  # result <- solve(p)
  # expect_equal(result$optimal_value, 0, tolerance = TOL)
  # expect_false(is.na(result$a))
})

test_that("Tests a problem with division", {
  obj <- Minimize(NormInf(A/5))
  p <- Problem(obj, list(A >= 5))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
})

