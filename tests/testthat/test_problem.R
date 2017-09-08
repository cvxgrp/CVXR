TOL <- 1e-6

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
  c1 <- matrix(rnorm(2), nrow = 1, ncol = 2)
  c2 <- matrix(rnorm(2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(c1 %*% x), list(x >= c2))
  # constants_ <- constants(p)
  # ref <- list(as.character(c1), as.character(c2))
  # mapply(function(c, r) { expect_equal(c, r) }, constants_, ref)
})

test_that("Test the size_metrics method", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  
  c1 <- matrix(rnorm(2), nrow = 2, ncol = 1)
  c2 <- matrix(rnorm(2), nrow = 1, ncol = 2)
  constants <- c(2, as.numeric(c2 %*% c1))
  
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + constants[1], c == constants[2]))
  
  # num_scalar variables
  n_variables <- size_metrics(p)@num_scalar_variables
  ref <- prod(size(a)) + prod(size(b)) + prod(size(c))
  expect_equal(n_variables, ref)
  
  # num_scalar_data
  n_data <- size_metrics(p)@num_scalar_data
  ref <- prod(size(p1)) + prod(size(p2)) + prod(size(p3)) + length(constants)  # 2 and c2 %*% c1 are both single scalar constants
  expect_equal(n_data, ref)
  
  # num_scalar_eq_constr
  n_eq_constr <- size_metrics(p)@num_scalar_eq_constr
  ref <- prod(dim(c2 %*% c1))
  expect_equal(n_eq_constr, ref)
  
  # num_scalar_leq_constr
  n_leq_constr <- size_metrics(p)@num_scalar_leq_constr
  ref <- prod(size(p3)) + prod(size(p2))
  expect_equal(n_leq_constr, ref)
  
  # max_data_dimension
  max_data_dim <- size_metrics(p)@max_data_dimension
  ref <- max(size(p3))
  expect_equal(max_data_dim, ref)
})

test_that("Test the solver_stats method", {
  prob <- Problem(Minimize(norm(x)), list(x == 0))
  result <- solve(prob, solver = "ECOS")
  # stats <- result$solver_stats
  # expect_true(stats$solve_time > 0)
  # expect_true(stats$setup_time > 0)
  # expect_true(stats$num_iters > 0)
  expect_true(result$solve_time > 0)
  expect_true(result$setup_time > 0)
  expect_true(result$num_iters > 0)
})

test_that("Test the get_problem_data method", {
  expect_error(get_problem_data(Problem(Maximize(Bool())), "ECOS"))
  
  data <- get_problem_data(Problem(Maximize(exp(a) + 2)), "SCS")
  dims <- data[["dims"]]
  expect_equal(dims[["ep"]], 1)
  expect_equal(length(data[["c"]]), 2)
  expect_equal(dim(data[["A"]]), c(3,2))
  
  data <- get_problem_data(Problem(Minimize(norm(x) + 3)), "ECOS")
  dims <- data[["dims"]]
  expect_equal(dims[["q"]], 3)
  expect_equal(length(data[["c"]]), 3)
  expect_equal(dim(data[["A"]]), c(0,3))
  expect_equal(dim(data[["G"]]), c(3,3))
})

test_that("Test silencing and enabling solver messages", {
  
})

test_that("Test registering other solve methods", {
  # Problem.register_solve("test", function(self) { 1 })
  # p <- Problem(Minimize(1))
  # result <- solve(p, method = "test")
  # expect_equal(result$optimal_value, 1)
  
  # test <- function(self, a, b = 2) { c(a, b) }
  # Problem.register_solve("test", test)
  # p <- Problem(Minimize(0))
  # result <- solve(p, 1, b = 3, method = "test")
  # expect_equal(result$optimal_value, c(1,3))
  # result <- solve(p, 1, method = "test")
  # expect_equal(result$optimal_value, c(1,2))
  # result <- solve(p, 1, method = "test", b = 4)
  # expect_equal(result$optimal_value, c(1,4))
})

test_that("Test that variables and constraints keep a consistent order", {
  num_solves <- 4
  vars_lists <- list()
  ineqs_lists <- list()
  var_ids_order_created <- list()
  for(k in 1:num_solves) {
    sum <- 0
    constraints <- list()
    var_ids <- c()
    for(i in 1:100) {
      var <- Variable(name = as.character(i))
      var_ids <- c(var_ids, var@id)
      sum <- sum + var
      constraints <- c(constraints, list(var >= i))
    }
    var_ids_order_created <- c(var_ids_order_created, list(var_ids))
    obj <- Minimize(sum)
    p <- Problem(obj, constraints)
    canon <- canonicalize(p)
    objective <- canon[[1]]
    constraints <- canon[[2]]
    sym_data <- SymData(objective, constraints, ECOS())
    
    # Sort by offset
    offsets <- sym_data@.var_offsets
    vars_ <- as.numeric(names(sort(offsets, decreasing = FALSE)))
    vars_lists <- c(vars_lists, list(vars_))
    ineqs_lists <- c(ineqs_lists, list(sym_data@.constr_map[[LEQ_MAP]]))
  }
  
  # Verify order of variables is consistent
  for(i in 1:num_solves)
    expect_equal(var_ids_order_created[[i]], vars_lists[[i]])
  
  for(i in 1:num_solves) {
    idx <- 1
    for(constr in ineqs_lists[[i]]) {
      var_tmp <- get_expr_vars(constr$expr)[[1]]
      var_id <- as.numeric(var_tmp[[1]])
      expect_equal(var_ids_order_created[[i]][idx], var_id)
      idx <- idx + 1
    }
  }
})

test_that("Test removing duplicate constraints objects", {
  eq <- (x == 2)
  le <- (x <= 2)
  obj <- 0
  
  test <- function(self) {
    canon <- canonicalize(self)
    objective <- canon[[1]]
    constraints <- canon[[2]]
    sym_data <- SymData(objective, constraints)   # TODO: What to use instead of CVXOPT?
    list(length(sym_data@constr_map[EQ_MAP]), length(sym_data@constr_map[LEQ_MAP]))
  }
  # Problem.register_solve("test", test)
  # p <- Problem(Minimize(obj), list(eq, eq, le, le))
  # result <- solve(p, method = "test")
  # expect_equal(result$optimal_value, c(1,1))
  
  # Internal constraints
  X <- Semidef(2)
  obj <- SumEntries(X + X)
  p <- Problem(Minimize(obj))
  # result <- solve(p, method = "test")
  # expect_equal(result$optimal_value, c(0,1))
  
  # Duplicates from non-linear constraints
  # exp <- norm(x, "2")
  # prob <- Problem(Minimize(0), list(exp <= 1, exp <= 2))
  # result <- solve(prob, method = "test")
  # expect_equal(result$optimal_value, c(0,4))
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
  value(p) <- 1
  
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
  # expect_equal(result$optimal_value, value(p))
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
  # expect_equal(result$optimal_value, value(p))
  # expect_equal(result$status, "INFEASIBLE")
  # expect_equal(result$optimal_value, Inf)
  # expect_true(result$optimal_value < 0)
  # expect_true(is.na(result$a))
  # expect_true(is.na(p@constraints[1]@dual_value))
  
  p <- Problem(Minimize(-a), list(a >= 2, a <= 1))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, value(p))
  # expect_equal(result$status, "INFEASIBLE")
  # expect_true(result$optimal_value, Inf)
  # expect_true(result$optimal_value > 0)
})

test_that("Test vector LP problems", {
  c <- Constant(matrix(c(1, 2), nrow = 2, ncol = 1))@value
  p <- Problem(Minimize(t(c) %*% x), list(x >= c))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 5, tolerance = TOL)
  # expect_equal(result$x, c(1, 2), tolerance = TOL)
  
  A <- Constant(rbind(c(3, 5), c(1, 2)))@value
  I <- Constant(rbind(c(1, 0), c(0, 1)))
  p <- Problem(Minimize(t(c) %*% x + a), list(A %*% x >= c(-1, 1), 4*I %*% z == x, z >= c(2,2), a >= 2))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 26, tolerance = 1e-3)
  # obj <- t(c) * result$x + result$a
  # expect_equal(obj[1,1], result$optimal_value, tolerance = TOL)
  # expect_equal(result$x, c(8,8), tolerance = 1e-3)
  # expect_equal(result$z, c(2,2), tolerance = 1e-3)
})

test_that("Test ECOS with no inequality constraints", {
  Tmat <- Constant(matrix(1, nrow = 2, ncol = 2))@value
  p <- Problem(Minimize(1), list(A == Tmat))
  # result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  # expect_equal(result$A, T, tolerance = TOL)
})

test_that("Test matrix LP problems", {
  Tmat <- Constant(matrix(1, nrow = 2, ncol = 2))@value
  p <- Problem(Minimize(1), list(A == Tmat))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  # expect_equal(result$A, T, tolerance = TOL)
  
  Tmat <- Constant(matrix(1, nrow = 2, ncol = 3)*2)@value
  c <- Constant(matrix(c(3,4), nrow = 2, ncol = 1))@value
  p <- Problem(Minimize(1), list(A >= Tmat %*% C, A == B, C == t(Tmat)))
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
  p <- Problem(Minimize(matrix(1, nrow = 1, ncol = 2) %*% (x + a + 1)), list(a + x >= c(1, 2)))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 5, tolerance = TOL)
})

test_that("Test parameter promotion", {
  a <- Parameter()
  value(a) <- 2
  exp <- cbind(c(1,2), c(3,4))*a
  expect_false(any(value(exp) - 2*cbind(c(1,2), c(3,4)) != 0))
})

test_that("test problems with parameters", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  p1@value <- 2
  p2@value <- -matrix(1, nrow = 3, ncol = 1)
  p3@value <- matrix(1, nrow = 4, ncol = 4)
  # result <- solve(p)
  # expect_equal(result$optimal_value, -6, tolerance = TOL)
  
  # p1@value <- NA
  # expect_error(solve(p))
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
  p <- Problem(Minimize(Norm2(5 + norm(z,1) + norm(x,1) + NormInf(x - z))), list(x >= c(2,3), z <= c(-1,-4), norm(x + z,2) <= 2))
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
    p <- Problem(Minimize(norm(x + z, 1)), list(x >= c(2,3), rbind(c(1,2), c(3,4)) %*% z == c(-1,4), norm(x + z, 2) <= 100))
    result <- solve(p, solver = solver)
    # expect_equal(result$optimal_value, 4, tolerance = acc)
    # expect_equal(result$x, c(4,3), tolerance = acc)
    # expect_equal(result$z, c(-4,1), tolerance = acc)
    
    # Dual values
    # expect_equal(result$constraints[1]$dual_value, c(0,1), tolerance = acc)
    # expect_equal(result$constraints[2]$dual_value, c(-1,0.5), tolerance = acc)
    # expect_equal(result$constraints[3]$dual_value, 0, tolerance = acc)
    
    Tmat <- matrix(1, nrow = 2, ncol = 3) * 2
    c <- matrix(c(3,4), nrow = 1, ncol = 2)
    p <- Problem(Minimize(1), list(A >= Tmat %*% C, A == B, C == t(Tmat)))
    result <- solve(p, solver = solver)
    
    # Dual values
    # expect_equal(result$constraints[1]$dual_value, rep(0,4), tolerance = acc)
    # expect_equal(result$constraints[2]$dual_value, rep(0,4), tolerance = acc)
    # expect_equal(result$constraints[3]$dual_value, rep(0,6), tolerance = acc)
  }
})

test_that("Test problems with indexing", {
  # Vector variables
  p <- Problem(Maximize(x[1,1]), list(x[1,1] <= 2, x[2,1] == 3))
  result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  # expect_equal(result$x, c(2,3))
  
  n <- 10
  Aloc <- matrix(0:(n^2-1), nrow = n, ncol = n)
  xloc <- Variable(n,n)
  p <- Problem(Minimize(SumEntries(xloc)), list(xloc == Aloc))
  result <- solve(p)
  answer <- n*n*(n*n+1)/2 - n*n
  # expect_equal(result$optimal_value, answer)
  
  # Matrix variables
  obj <- A[1,1] + A[1,2] + A[2,2] + A[2,1]
  p <- Problem(Maximize(obj), list(A <= rbind(c(1,-2), c(-3,4))))
  result <- solve(p)
  # expect_equal(result$optimal_value, 0, tolerance = TOL)
  # expect_equal(result$A, c(1,-2,-3,4))
  
  # Indexing arithmetic expressions
  exp <- rbind(c(1,2), c(3,4)) %*% z + x
  p <- Problem(Minimize(exp[2,1]), list(x == z, z == c(1,2)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  # expect_equal(result$x, result$z, tolerance = TOL)
})

test_that("Test problems with slicing", {
  p <- Problem(Maximize(SumEntries(C)), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(result$C, 2*c(1,2,2))
  
  p <- Problem(Maximize(SumEntries(C[seq(1,3,2),2])), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 3, tolerance = TOL)
  # expect_equal(result$C[seq(1,3,2),2], c(1,2))
  
  p <- Problem(Maximize(SumEntries((C[1:2,] + A)[,1:2])), list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3,
                                                             (A + B)[,2] == 2, B == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  # expect_equal(result$C[1:2,], c(1,2,1,2), tolerance = TOL)
  # expect_equal(result$A, c(2,2,1,1), tolerance = TOL)
  
  p <- Problem(Maximize(matrix(c(3,4), nrow = 1, ncol = 2) %*% (C[1:2,] + A)[,1]),
               list(C[2:3,] <= 2, C[1,] == 1, matrix(c(1,2), nrow = 1, ncol = 2) %*% (A + B)[,1] == 3,
                    (A + B)[,2] == 2, B == 1, 3*A[,1] <= 3))
  result <- solve(p)
  # expect_equal(result$optimal_value, 12, tolerance = TOL)
  # expect_equal(result$C[1:2,1], c(1,2), tolerance = TOL)
  # expect_equal(result$A, c(1,-0.5,1,1), tolerance = TOL)
  
  p <- Problem(Minimize(Norm2(C[1:2,] + A)[,1]), list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3, (A + B)[,2] == 2, B == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 3, tolerance = TOL)
  # expect_equal(result$C[1:2,1], c(1,-2), tolerance = TOL)
  # expect_equal(result$A, c(2,2,1,1), tolerance = TOL)
  
  # Transpose of slice
  p <- Problem(Maximize(SumEntries(C)), list(t(C[2:3,]) <= 2, t(C[1,]) == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(result$C, 2*c(1,2,2))
})

test_that("Test the VStack atom", {
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% VStack(x, y)), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 15, tolerance = TOL)
  
  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% VStack(x, x)), list(x == c(1,2)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  
  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(SumEntries(VStack(A, C))), list(A >= 2*c, C == -2))
  result <- solve(p)
  # expect_equal(result$optimal_value, -4, tolerance = TOL)
  
  c <- matrix(1, nrow = 1, ncol = 2)
  p <- Problem(Minimize(SumEntries(VStack(c %*% A, c %*% B))), list(A >= 2, B == -2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 0, tolerance = TOL)
  
  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% VStack(Square(a), sqrt(b))), list(a == 2, b == 16))
  expect_error(solve(p))
})

test_that("Test the HStack atom", {
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% t(HStack(t(x), t(y)))), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 15, tolerance = TOL)
  
  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% t(HStack(t(x), t(x)))), list(x == c(1,2)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  
  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(SumEntries(HStack(t(A), t(C)))), list(A >= 2*c, C == -2))
  result <- solve(p)
  # expect_equal(result$optimal_value, -4, tolerance = TOL)
  
  D <- Variable(3,3)
  expr <- HStack(C, D)
  p <- Problem(Minimize(expr[1,2] + SumEntries(HStack(expr, expr))), list(C >= 0, D >= 0, D[1,1] == 2, C[1,2] == 3))
  result <- solve(p)
  # expect_equal(result$optimal_value, 13, tolerance = TOL)
  
  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% t(HStack(t(Square(a)), t(sqrt(b))))), list(a == 2, b == 16))
  expect_error(solve(p))
})

test_that("Test using a CVXR expression as an objective", {
  expect_error(Problem(x+2))
})

test_that("Test variable transpose", {
  p <- Problem(Minimize(SumEntries(x)), list(t(x) >= matrix(c(1,2), nrow = 1, ncol = 2)))
  result <- solve(p)
  # expect_equal(result$optimal_value, 3, tolerance = TOL)
  # expect_equal(result$x, c(1,2))
  
  p <- Problem(Minimize(SumEntries(C)), list(matrix(c(1,1), nrow = 1, ncol = 2) %*% t(C) >= matrix(0:2, nrow = 1, ncol = 3)))
  result <- solve(p)
  # value <- result$C
  
  constraints <- lapply(1:3, function(i) { 1*C[i,1] + 1*C[i,2] >= i })
  p <- Problem(Minimize(SumEntries(C)), constraints)
  result2 <- solve(p)
  # expect_equal(result$optimal_value, result2$optimal_value)
  # expect_equal(result$C, value)
  
  p <- Problem(Minimize(A[1,2] - t(A)[2,1]), list(A == cbind(c(1,2), c(3,4))))
  result <- solve(p)
  # expect_equal(result$optimal_value, 0, tolerance = TOL)
  
  exp <- t(-x)
  p <- Problem(Minimize(SumEntries(x)), list(t(-x) <= 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, -2, tolerance = TOL)
  
  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(MaxElemwise(t(c), 2, 2 + t(c))[2]))
  result <- solve(p)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(SumEntries(t(MaxElemwise(c, 2, 2+c))[,1])))
  result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  
  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(SumEntries(t(Square(t(c)))[,1])))
  result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  
  # Slice of transpose
  p <- Problem(Maximize(SumEntries(C)), list(t(C)[,2:3] <= 2, t(C)[,1] == 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(result$C, 2*c(1,2,2))
})

test_that("Test multiplication on the left by a non-constant", {
  c <- matrix(c(1,2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% A %*% c), list(A >= 2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 18, tolerance = TOL)
  
  p <- Problem(Minimize(a*2), list(a >= 2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 4, tolerance = TOL)
  
  p <- Problem(Minimize(t(x) %*% c), list(x >= 2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 6, tolerance = TOL)
  
  p <- Problem(Minimize((t(x) + t(z)) %*% c), list(x >= 2, z >= 1))
  result <- solve(p)
  # expect_equal(result$optimal_value, 9, tolerance = TOL)
})

test_that("Test redundant constraints", {
  obj <- Minimize(SumEntries(x))
  constraints <- list(x == 2, x == 2, t(x) == 2, x[1] == 2)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, 4, tolerance = TOL)
  
  obj <- Minimize(SumEntries(Square(x)))
  constraints <- list(x == x)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "ECOS")
  # expect_equal(result$optimal_value, 0, tolerance = TOL)
})

test_that("Test that symmetry is enforced", {
  p <- Problem(Minimize(LambdaMax(A)), list(A >= 2))
  result <- solve(p)
  # expect_equal(result$A, t(result$A), tolerance = 1e-3)
  
  p <- Problem(Minimize(LambdaMax(A)), list(A == cbind(c(1,2), c(3,4))))
  result <- solve(p)
  expect_equal(tolower(result$status), "infeasible")
})

test_that("Test SDP", {
  # Ensure SDP constraints enforce transpose
  obj <- Maximize(A[2,1] - A[1,2])
  p <- Problem(obj, list(LambdaMax(A) <= 100, A[1,1] == 2, A[2,2] == 2, A[2,1] == 2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 0, tolerance = 1e-3)
})

test_that("Test getting values for expressions", {
  diff_exp <- x - z
  inf_exp <- NormInf(diff_exp)
  sum_entries_exp <- 5 + norm(z,1) + norm(x,1) + inf_exp
  constr_exp <- Norm2(x + z)
  obj <- Norm2(sum_entries_exp)
  p <- Problem(Minimize(obj), list(x >= c(2,3), z <= c(-1,-4), constr_exp <= 2))
  result <- solve(p)
  # expect_equal(result$optimal_value, 22, tolerance = TOL)
  # expect_equal(result$x, c(2,3))
  # expect_equal(result$z, c(-1,-4))
  
  # Expression values
  # expect_equal(value(diff_exp, result), result$x - result$z, tolerance = TOL)
  # expect_equal(value(inf_exp, result), norm(result$x - result$z, "I"))
  # expect_equal(value(sum_entries_exp, result), 5 + norm(result$z, "1") + norm(result$x, "1") + norm(result$x - result$z, "I"))
  # expect_equal(value(constr_exp, result), norm(result$x + result$z, "2"))
  # expect_equal(value(obj, result), result$optimal_value)
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

test_that("Tests problems with MulElemwise", {
  c <- cbind(c(1,-1), c(2,-2))
  expr <- MulElemwise(c, A)
  obj <- Minimize(NormInf(expr))
  p <- Problem(obj, list(A == 5))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(value(expr, result), c(5,-5) + c(10,-10), tolerance = TOL)
  
  # Test with a sparse matrix
  c <- sparseMatrix(i = c(1,2), j = c(1,1), x = c(1,2))
  expr <- MulElemwise(c, x)
  obj <- Minimize(NormInf(expr))
  p <- Problem(obj, list(x == 5))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(value(expr, result), c(5,10))
  
  # Test promotion
  c <- cbind(c(1,-1), c(2,-2))
  expr <- MulElemwise(c, a)
  obj <- Minimize(NormInf(expr))
  p <- Problem(obj, list(a == 5))
  result <- solve(p)
  # expect_equal(result$optimal_value, 10, tolerance = TOL)
  # expect_equal(value(expr, result), c(5,-5) + c(10,-10))
})

test_that("Tests that errors occur when you use an invalid solver", {
  expect_error(solve(Problem(Minimize(Bool())), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(LambdaMax(a))), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(a)), solver = "SCS"))
})

test_that("Tests problems with Reshape", {
  # Test on scalars
  expect_equal(value(Reshape(1,1,1)), 1)
  
  # Test vector to matrix
  x <- Variable(4)
  mat <- cbind(c(1,-1), c(2,-2))
  vec <- matrix(1:4)
  vec_mat <- cbind(c(1,2), c(3,4))
  expr <- Reshape(x,2,2)
  obj <- Minimize(SumEntries(mat %*% expr))
  prob <- Problem(obj, list(x == vec))
  result <- solve(prob)
  # expect_equal(result$optimal_value, sum(mat %*% vec_mat))
  
  # Test on matrix to vector
  c <- 1:4
  expr <- Reshape(A,4,1)
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == rbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  # expect_equal(result$optimal_value, 20, tolerance = TOL)
  # expect_equal(value(expr, result), c(-1,-2,3,4))
  # expect_equal(value(Reshape(expr,2,2), result), c(-1,-2,3,4))
  
  # Test matrix to matrix
  expr <- Reshape(C,2,3)
  mat <- rbind(c(1,-1), c(2,-2))
  C_mat <- rbind(c(1,4), c(2,5), c(3,6))
  obj <- Minimize(SumEntries(mat %*% expr))
  prob <- Problem(obj, list(C == C_mat))
  result <- solve(prob)
  reshaped = matrix(C_mat, nrow = 2, ncol = 3, byrow = FALSE)
  # expect_equal(result$optimal_value, sum(mat %*% reshaped), tolerance = TOL)
  # expect_equal(value(expr, result), C_mat, tolerance = TOL)
  
  # Test promoted expressions
  c <- cbind(c(1,-1), c(2,-2))
  expr <- Reshape(c * a,1,4)
  obj <- Minimize(expr %*% (1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob)
  # expect_equal(result$optimal_value, -6, tolerance = TOL)
  # expect_equal(value(expr, result), 2*c, tolerance = TOL)
  
  expr <- Reshape(c * a,4,1)
  obj <- Minimize(t(expr) %*% (1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob)
  # expect_equal(result$optimal_value, -6, tolerance = TOL)
  # expect_equal(value(expr, result), 2*c, tolerance = TOL)
})

test_that("Tests problems with Vec", {
  c <- 1:4
  expr <- Vec(A)
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == cbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  # expect_equal(result$optimal_value, 20, tolerance = TOL)
  # expect_equal(value(expr, result), c(-1,-2,3,4))
})

test_that("Test a problem with Diag", {
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(Diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == Semidef(3))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  # expect_equal(result$optimal_value, 0.583151, tolerance = 1e-2)
})

test_that("Test that the presolver removes constraints with no variables", {
  x <- Variable()
  obj <- Maximize(sqrt(x))
  prob <- Problem(obj, list(Constant(2) <= 2))
  data <- get_problem_data(prob, "ECOS")
  A <- data[["A"]]
  G <- data[["G"]]
  if(nrow(A) > 0)
    tmp <- apply(A, 1, function(row) { expect_true(sum(row != 0) > 0) })
  if(nrow(G) > 0)
    tmp <- apply(G, 1, function(row) { expect_true(sum(row != 0) > 0) })
})

test_that("Test presolve parameters", {
  # # Test with parameters
  # gamma <- Parameter(sign = "positive")
  # x <- Variable()
  # obj <- Minimize(x)
  # value(gamma) <- 0
  # prob <- Problem(obj, list(gamma == 1, x >= 0))
  # result <- solve(prob, solver = "SCS")
  # expect_equal(tolower(result$status), "infeasible")
  # 
  # value(gamma) <- 1
  # prob <- Problem(obj, list(gamma == 1, x >= 0))
  # result <- solve(prob, solver = "SCS")
  # expect_equal(tolower(result$status), "optimal")
})

test_that("Test that expressions with parameters are updated properly", {
  # x <- Variable()
  # y <- Variable()
  # x0 <- Parameter()
  # 
  # # Initial guess for x
  # value(x0) <- 2
  # 
  # # Make the constraints x^2 - y == 0
  # xSquared <- x0*x0 + 2*x0*(x-x0)
  # g <- xSquared - y
  # 
  # # Set up the problem
  # obj <- abs(x - 1)
  # prob <- Problem(Minimize(obj), list(g == 0))
  # result <- solve(prob)
  # 
  # value(x0) <- 1
  # xSquared <- x0*x0 + 2*x0*(x-x0)
  # g <- xSquared - y
  # prob <- Problem(Minimize(obj), list(g == 0))
  # result <- solve(prob)
  # # expect_equal(result$g, 0, tolerance = TOL)
  # 
  # # Test multiplication
  # prob <- Problem(Minimize(x0*x), list(x == 1))
  # value(x0) <- 2
  # result <- solve(prob)
  # x0@value <- 1
  # result <- solve(prob)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
})

test_that("Test interaction of caching with changing constraints", {
  prob <- Problem(Minimize(a), list(a == 2, a >= 1))
  result <- solve(prob)
  # expect_equal(result$optimal_value, 2, tolerance = TOL)
  
  prob@constraints[[1]] = (a == 1)
  result <- solve(prob)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
})

test_that("Test positive definite constraints", {
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(Diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == t(C), C %>>% 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  # expect_equal(result$optimal_value, 0.583151, tolerance = 1e-2)
  
  C <- Variable(2,2)
  obj <- Maximize(C[1,2])
  constraints <- list(C == 1, C %>>% rbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(tolower(result$status), "infeasible")
  
  C <- Symmetric(2,"2")
  obj <- Minimize(C[1,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(tolower(result$status), "unbounded")
})

test_that("Test the duals of PSD constraints", {
  C <- Symmetric(2,"2")
  obj <- Maximize(C[1,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  # expect_equal(result$optimal_value, 2, tolerance = 1e-4)
  
  psd_constr_dual <- dual_value(constraints[[1]])
  C <- Symmetric(2,"2")
  X <- Semidef(2)
  obj <- Maximize(C[1,1])
  constraints <- list(X == cbind(c(2,0), c(0,2)) - C)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  # expect_equal(constraints[1]@dual_value, psd_constr_dual)
  
  # Test dual values with SCS that have off-diagonal entries
  C <- Symmetric(2,"2")
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)), C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  # expect_equal(result$optimal_value, 4, tolerance = 1e-3)
  
  psd_constr_dual <- dual_value(constraints[[1]])
  C <- Symmetric(2,"2")
  X <- Semidef(2)
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(X == cbind(c(2,0), c(0,2)) - C, C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  # expect_equal(dual_value(constraints[[1]]), psd_constr_dual, tolerance = 1e-3)
})

test_that("Test GeoMean", {
  require(gmp)
  x <- Variable(2)
  cost <- GeoMean(x)
  prob <- Problem(Maximize(cost), list(x <= 1))
  result <- solve(prob)
  # expect_equal(result$optimal_value, 1, tolerance = TOL)
  
  prob <- Problem(Maximize(cost), list(sum(x) <= 1))
  result <- solve(prob)
  # expect_equal(result$x, c(0.5,0.5))
  
  x <- Variable(3,3)
  expect_error(GeoMean(x))
  
  x <- Variable(3,1)
  g <- GeoMean(x)
  # expect_equal(g@w, as.bigq(1,3)*3, tolerance = TOL)
  
  x <- Variable(1,5)
  g <- GeoMean(x)
  # expect_equal(g@w, as.bigq(1,5)*5, tolerance = TOL)
  
  # Check that we get the right answer for max GeoMean(x) s.t. sum(x) <= 1
  p <- c(0.07, 0.12, 0.23, 0.19, 0.39)
  
  short_geo_mean <- function(x, p) {
    p <- as.numeric(p)/sum(p)
    x <- as.numeric(x)
    prod(x^p)
  }
  
  x <- Variable(5)
  prob <- Problem(Maximize(GeoMean(x, p)), list(sum(x) <= 1))
  result <- solve(prob)
  # x <- as.numeric(result$x)
  # x_true <- p/sum(p)
  
  # expect_equal(result$optimal_value, value(GeoMean(list(x), p)), tolerance = TOL)
  # expect_equal(result$optimal_value, short_geo_mean(x, p), tolerance = TOL)
  # expect_equal(x, x_true, tolerance = 1e-3)
  
  # Check that we get the right answer for max GeoMean(x) s.t. norm(x) <= 1
  x <- Variable(5)
  prob <- Problem(Maximize(GeoMean(x, p)), list(norm(x) <= 1))
  result <- solve(prob)
  # x <- as.numeric(result$x)
  # x_true <- sqrt(p/sum(p))
  
  # expect_true(result$optimal_value, value(GeoMean(list(x), p)), tolerance = TOL)
  # expect_true(result$optimal_value, short_geo_mean(x, p), tolerance = TOL)
  # expect_true(x, x_true, tolerance = 1e-3)
  
  # The following 3 tests check VStack and HStack input to GeoMean
  # The following 3 formulations should be equivalent
  n <- 5
  x_true <- rep(1,n)
  x <- Variable(n)
  
  result <- solve(Problem(Maximize(GeoMean(x)), list(x <= 1)))
  # xval <- as.numeric(result$x)
  # expect_equal(xval, x_true, tolerance = 1e-3)
  
  args <- list()
  for(i in 1:n) 
    args <- c(args, x[i])
  
  y <- do.call(VStack, args)
  result <- solve(Problem(Maximize(GeoMean(y)), list(x <= 1)))
  # xval <- as.vector(result$x)
  # expect_equal(xval, x_true, tolerance = 1e-3)
  
  y <- do.call(HStack, args)
  result <- solve(Problem(Maximize(GeoMean(y)), list(x <= 1)))
  # xval <- as.vector(result$x)
  # expect_equal(xval, x_true, tolerance = 1e-3)
})

test_that("Test Pnorm", {
  x <- Variable(3, name = "x")
  avec <- c(1.0, 2, 3)
  
  # TODO: Add -1, 0.5, 0.3, -2.3 and testing positivity constraints
  
  for(p in c(1, 1.6, 1.3, 2, 1.99, 3, 3.7, Inf)) {
    prob <- Problem(Minimize(Pnorm(x, p = p)), list(t(x) %*% avec >= 1))
    result <- solve(prob)
    
    # Formula is true for any a >= 0 with p > 1
    if(p == Inf)
      x_true <- rep(1, length(avec))/sum(avec)
    else if(p == 1) {
      # Only works for the particular a = c(1,2,3)
      x_true <- c(0,0,1.0/3)
    } else
      x_true <- avec^(1.0/(p-1)) / as.numeric(avec %*% (avec^(1.0/(p-1))))
    
    x_alg <- as.vector(result$x)
    # expect_equal(x_alg, x_true, tolerance = 1e-3)
    # expect_equal(result$optimal_value, norm(x_alg, p))
    # expect_equal(norm(x_alg, p), value(Pnorm(x_alg, p), result))
  }
})

test_that("Test Pnorm concave", {
  x <- Variable(3, name = "x")
  
  # Test positivity constraints
  a <- c(-1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(SumEntries(abs(x-a))), list(Pnorm(x,p) >= 0))
    result <- solve(prob)
    
    # expect_equal(result$optimal_value, 1)
  }
  
  a <- c(1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(SumEntries(abs(x-a))), list(Pnorm(x,p) >= 0))
    result <- solve(prob)
    
    # expect_equal(result$optimal_value, 0)
  }
})

test_that("Test Power", {
  x <- Variable()
  prob <- Problem(Minimize(Power(x, 1.7) + Power(x, -2.3) - Power(x, 0.45)))
  result <- solve(prob)
  # x <- result$x
  # expect_true(abs(1.7*x^0.7 - 2.3*x^-3.3 - 0.45*x^-0.55) <= 1e-3)
})

test_that("Test a problem with MulElemwise by a scalar", {
  Tnum <- 10
  Jnum <- 20
  rvec <- matrix(rnorm(Tnum*Jnum), nrow = Tnum, ncol = Jnum)
  dy <- matrix(rnorm(2*Tnum), nrow = 2*Tnum, ncol = 1)
  theta <- Variable(Jnum)
  
  delta <- 1e-3
  loglambda <- rvec %*% theta   # rvec: TxJ regressor matrix, theta: (Jx1) cvx variable
  a <- MulElemwise(dy[1:Tnum], loglambda)  # size (Tx1)
  b1 <- exp(loglambda)
  b2 <- MulElemwise(delta, b1)
  cost <- -a + b1
  
  cost <- -a + b2  # size (Tx1)
  prob <- Problem(Minimize(SumEntries(cost)))
  result <- solve(prob, solver = "SCS")
  
  # obj <- Minimize(SumEntries(MulElemwise(2, result$x)))
  # prob <- Problem(obj, list(result$x == 2))
  # result <- solve(prob)
  # expect_equal(result$optimal_value, 8, tolerance = TOL)
})
