context("test-g01-dgp2dcp")

Dgp2Dcp <- function(problem) { new("Dgp2Dcp", problem = problem) }

test_that("test unconstrained monomial", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod <- x*y
  dgp <- Problem(Minimize(prod))
  dgp2dcp <- Dgp2Dcp(dgp)
  
  dcp <- reduce(dgp2dcp)
  expect_equal(class(dcp@objective@expr), "AddExpression")
  expect_equal(length(dcp@objective@expr@args), 2)
  expect_equal(class(dcp@objective@expr@args[[1]]), "Variable")
  expect_equal(class(dcp@objective@expr@args[[2]]), "Variable")
  opt <- solve(dcp)
  
  # dcp is solved in log-space, so it is unbounded below
  # (since the OPT for dgp is 0 + epsilon).
  expect_equal(opt$value, -Inf)
  expect_equal(opt$status, "unbounded")
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
  expect_equal(dgp@value, 0.0)
  expect_equal(dgp@status, "unbounded")
  dgp <- clear_solution(dgp)
  opt <- solve(dgp, gp = TRUE)
  expect_equal(opt$value, 0.0)
  expect_equal(opt$status, "unbounded")
  
  dgp <- Problem(Maximize(prod))
  dgp2dcp <- Dgp2Dcp(dgp)
  dcp <- reduce(dgp2dcp)
  opt <- solve(dcp)
  expect_equal(opt$value, Inf)
  expect_equal(opt$status, "unbounded")
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
  expect_equal(dgp@value, Inf)
  expect_equal(dgp@status, "unbounded")
  dgp <- clear_solution(dgp)
  opt <- solve(dgp, gp = TRUE)
  expect_equal(opt$value, Inf)
  expect_equal(opt$status, "unbounded")
})

test_that("test basic equality constraint", {
  x <- Variable(pos = TRUE)
  dgp <- Problem(Minimize(x), list(x == 1.0))
  dgp2dcp <- Dgp2Dcp(dgp)
  
  dcp <- reduce(dgp2dcp)
  expect_equal(class(dcp@objective@expr), "Variable")
  opt <- solve(dcp)
  expect_equal(opt$value, 0.0, tolerance = TOL)
  expect_equal(value(variables(dcp)[[1]]), 0.0, tolerance = TOL)
  
  dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
  expect_equal(dgp@value, 1.0, tolerance = TOL)
  expect_equal(value(x), 1.0, tolerance = TOL)
  dgp <- clear_solution(dgp)
  result <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
})

test_that("test basic GP", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  constraints <- list(2*x*y + 2*x*z + 2*y*z <= 1.0, x >= 2*y)
  problem <- Problem(Minimize(1/(x*y*z)), constraints)
  result <- solve(problem, gp = TRUE)
  expect_equal(15.59, result$value, tolerance = 1e-2)
})

test_that("test max_elemwise", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  prod1 <- x*y^(0.5)
  prod2 <- 3.0*x*y^(0.5)
  obj <- Minimize(max_elemwise(prod1, prod2))
  constr <- list(x == 1.0, y == 4.0)
  
  dgp <- Problem(obj, constr)
  dgp2dcp <- Dgp2Dcp(dgp)
  dcp <- reduce(dgp2dcp)
  opt <- solve(dcp)
  unpack(dgp, retrieve(dgp2dcp, opt$solution))
  expect_equal(dgp@value, 6.0, tolerance = TOL)
  expect_equal(value(x), 1.0, tolerance = TOL)
  expect_equal(value(y), 4.0, tolerance = TOL)
  dgp <- clear_solution(dgp)
  opt <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 6.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
})

# TODO: Test keepdims.
test_that("test prod_entries", {
  X <- matrix(0:11, nrow = 4, ncol = 3)
  expect_equal(prod(X), value(prod_entries(X)))
  expect_equal(apply(X, 1, prod), value(prod_entries(X, axis = 1)))
  expect_equal(apply(X, 2, prod), value(prod_entries(X, axis = 2)))
  
  prod <- prod_entries(X)
  X_canon <- Dgp2Dcp.prod_canon(prod, prod@args)[[1]]
  expect_equal(sum(X), value(X_canon))
  
  prod <- prod_entries(X, axis = 1)
  X_canon <- Dgp2Dcp.prod_canon(prod, prod@args)[[1]]
  expect_equal(apply(X, sum, axis = 1), value(X_canon))
  
  prod <- prod_entries(X, axis = 2)
  X_canon <- Dgp2Dcp.prod_canon(prod, prod@args)[[1]]
  expect_equal(apply(X, sum, axis = 2), value(X_canon))
  
  X <- matrix(0:11, nrow = 12, ncol = 1)
  expect_equal(prod(X), value(prod_entries(X)))
  
  prod <- prod_entries(X)
  X_canon <- Dgp2Dcp.prod_canon(prod, prod@args)[[1]]
  expect_equal(sum(X), value(X_canon))
  
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  posy1 <- x*y^(0.5) + 3.0*x*y^(0.5)
  posy2 <- x*y^(0.5) + 3.0*x^2*y^(0.5)
  expect_true(is_log_log_convex(prod(posy1, posy2)))
  expect_false(is_log_log_concave(prod(posy1, posy2)))
  expect_false(is_dgp(prod(posy1, 1/posy1)))
  
  m <- x*y^(0.5)
  expect_true(is_log_log_affine(prod(m, m)))
  expect_true(is_log_log_concave(prod(m, 1/posy1)))
  expect_false(is_log_log_convex(prod(m, 1/posy1)))
})

test_that("test max_entries", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  prod1 <- x*y^(0.5)
  prod2 <- 3.0*x*y^(0.5)
  obj <- Minimize(max_entries(hstack(list(prod1, prod2))))
  constr <- list(x == 1.0, y == 4.0)
  
  dgp <- Problem(obj, constr)
  dgp2dcp <- Dgp2Dcp(dgp)
  dcp <- reduce(dgp2dcp)
  opt <- solve(dcp)
  dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
  expect_equal(dgp@value, 6.0)
  expect_equal(value(x), 1.0)
  expect_equal(value(y), 4.0)
  dgp <- clear_solution(dgp)
  result <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 6.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
})

test_that("test min_elemwise", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  prod1 <- x*y^(0.5)
  prod2 <- 3.0*x*y^(0.5)
  posy <- prod1 + prod2
  obj <- Maximize(min_elemwise(prod1, prod2, 1/posy))
  constr <- list(x == 1.0, y == 4.0)
  
  dgp <- Problem(obj, constr)
  result <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 1.0/(2.0 + 6.0), tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
})

test_that("test min_entries", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  prod1 <- x*y^(0.5)
  prod2 <- 3.0*x*y^(0.5)
  posy <- prod1 + prod2
  obj <- Maximize(min_entries(hstack(list(prod1, prod2, 1/posy))))
  constr <- list(x == 1.0, y == 4.0)
  
  dgp <- Problem(obj, constr)
  result <- solve(dgp, gp = TRUE)
  expect_equal(result$value, 1.0/(2.0 + 6.0), tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
})

# TODO: Enable test once sum_largest is implemented.
# test_that("test sum_largest", {
#   x <- Variable(4, pos = TRUE)
#   obj <- Minimize(sum_largest(x, 3))
#   constr <- list(x[1]*x[2]*x[3]*x[4] >= 16)
#   dgp <- Problem(obj, constr)
#   dgp2dcp <- Dgp2Dcp(dgp)
#   dcp <- reduce(dgp2dcp)
#   result <- solve(dcp)
#   dgp <- unpack(dgp, retrieve(dgp2dcp, result$solution))
#   opt <- 6.0
#   expect_equal(dgp@value, opt, tolerance = TOL)
#   expect_equal(value(x[1]*x[2]*x[3]*x[4]), 16, tolerance = 1e-2)
#   dgp <- clear_solution(dgp)
#   result <- solve(dgp, gp = TRUE)
#   expect_equal(result$value, opt, tolerance = TOL)
#   expect_equal(value(x[1]*x[2]*x[3]*x[4]), 16, tolerance = 1e-2)
#   
#   # An unbounded problem.
#   x <- Variable(4, pos = TRUE)
#   y <- Variable(pos = TRUE)
#   obj <- Minimize(sum_largest(x, 3)*y)
#   constr <- list(x[1]*x[2]*x[3]*x[4] >= 16)
#   dgp <- Problem(obj, constr)
#   dgp2dcp <- Dgp2Dcp(dgp)
#   dcp <- reduce(dgp2dcp)
#   opt <- solve(dcp)
#   expect_equal(dcp@value, -Inf)
#   dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
#   expect_equal(dgp@value, 0.0, tolerance = TOL)
#   expect_equal(dgp@status, "unbounded")
#   dgp <- clear_solution(dgp)
#   result <- solve(dgp, gp = TRUE)
#   expect_equal(result$value, 0.0, tolerance = TOL)
#   expect_equal(result$status, "unbounded")
#   
#   # Another unbounded problem.
#   x <- Variable(2, pos = TRUE)
#   obj <- Minimize(sum_largest(x, 1))
#   dgp <- Problem(obj)
#   dgp2dcp <- Dgp2Dcp(dgp)
#   dcp <- reduce(dgp2dcp)
#   opt <- solve(dcp)
#   expect_equal(result$value, -Inf)
#   dgp <- unpack(dgp, retrieve(dgp2dcp, opt$solution))
#   expect_equal(dgp@value, 0.0, tolerance = TOL)
#   expect_equal(dgp@status, "unbounded")
#   dgp <- clear_solution(dgp)
#   result <- solve(dgp, gp = TRUE)
#   expect_equal(dgp@value, 0.0, tolerance = TOL)
#   expect_equal(dgp@status, "unbounded")
#   
#   # Composition with posynomials.
#   x <- Variable(4, pos = TRUE)
#   obj <- Minimize(sum_largest(hstack(list(3*x[1]^(0.5) * x[2]^(0.5), x[1]*x[2] + 0.5*x[2]*x[4]^3, x[3])), 2))
#   constr <- list(x[1]*x[2] >= 16)
#   dgp <- Problem(obj, constr)
#   dgp2dcp <- Dgp2Dcp(dgp)
#   dcp <- reduce(dgp2dcp)
#   result <- solve(dcp)
#   dgp <- unpack(dgp, retrieve(dgp2dcp, result$solution))
#   
#   # opt = 3 * sqrt(4) * sqrt(4) + (4 * 4 + 0.5 * 4 * epsilon) = 28
#   opt <- 28.0
#   expect_equal(dgp@value, opt, tolerance = 1e-2)
#   expect_equal(value(x[1]*x[2]), 16.0, tolerance = 1e-2)
#   expect_equal(value(x[4]), 0.0, tolerance = 1e-2)
#   dgp <- clear_solution(dgp)
#   result <- solve(dgp, gp = TRUE)
#   expect_equal(result$value, opt, tolerance = 1e-2)
#   expect_equal(value(x[1]*x[2]), 16.0, tolerance = 1e-2)
#   expect_equal(value(x[4]), 0.0, tolerance = 1e-2)
# })

test_that("test div", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  p <- Problem(Minimize(x*y), list(y/3 <= x, y >= 1))
  result <- solve(p, gp = TRUE)
  expect_equal(result$value, 1.0/3.0, tolerance = TOL)
  expect_equal(result$getValue(y), 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0/3.0, tolerance = TOL)
})

test_that("test geo_mean", {
  x <- Variable(3, pos = TRUE)
  p <- c(1, 2, 0.5)
  expr <- geo_mean(x, p)
  dgp <- Problem(Minimize(expr))
  dgp2dcp <- Dgp2Dcp(dgp)
  dcp <- reduce(dgp2dcp)
  result <- solve(dcp)
  expect_equal(result$value, -Inf)
  dgp <- unpack(dgp, retrieve(dgp2dcp, result$solution))
  expect_equal(dgp@value, 0.0)
  expect_equal(dgp@status, "unbounded")
  dgp <- clear_solution(dgp)
  result <- solve(dgp, gp = TRUE)
  expect_equal(dgp@value, 0.0)
  expect_equal(dgp@status, "unbounded")
})

test_that("test solving non-dgp problem raises error", {
  problem <- Problem(Minimize(-1.0*Variable()))
  expect_error(solve(problem, gp = TRUE))
  result <- solve(problem)
  expect_equal(result$status, "unbounded")
  expect_equal(result$value, -Inf)
})

test_that("test solving non-dcp problem raises error", {
  problem <- Problem(Minimize(Variable(pos = TRUE) * Variable(pos = TRUE)))
  expect_error(solve(problem))
  result <- solve(problem, gp = TRUE)
  expect_equal(result$status, "unbounded", tolerance = TOL)
  expect_equal(result$value, 0.0, tolerance = TOL)
})

test_that("test add_canon", {
  X <- Constant(rbind(1:3, 4:6))
  Y <- Constant(rbind(2:4, 5:7))
  Z <- X + Y
  canon <- Dgp2Dcp.add_canon(Z, Z@args)
  canon_matrix <- canon[[1]]
  constraints <- canon[[2]]
  expect_equal(length(constraints), 0)
  expect_equal(dim(canon_matrix), dim(Z))
  expected <- log(exp(value(X)) + exp(value(Y)))
  expect_equal(expected, value(canon_matrix), tolerance = TOL)
  
  # Test promotion.
  X <- Constant(rbind(1:3, 4:6))
  y <- Constant(2.0)
  Z <- X + y
  canon <- Dgp2Dcp.add_canon(Z, Z@args)
  canon_matrix <- canon[[1]]
  constraints <- canon[[2]]
  expect_equal(length(constraints), 0)
  expect_equal(dim(canon_matrix), dim(Z))
  expected <- log(exp(value(X)) + exp(value(y)))
  expect_equal(expected, value(canon_matrix), tolerance = TOL)
})

test_that("test matmul_canon", {
  X <- Constant(rbind(1:3, 4:6))
  Y <- Constant(matrix(1:3, nrow = 3))
  Z <- matmul(X, Y)
  canon <- Dgp2Dcp.mulexpression_canon(Z, Z@args)
  canon_matrix <- canon[[1]]
  constraints <- canon[[2]]
  expect_equal(length(constraints), 0)
  expect_equal(dim(canon_matrix), c(2,1))
  first_entry <- log(exp(2.0) + exp(4.0) + exp(6.0))
  second_entry <- log(exp(5.0) + exp(7.0) + exp(9.0))
  expect_equal(first_entry, value(canon_matrix[1,1]), tolerance = TOL)
  expect_equal(second_entry, value(canon_matrix[2,1]), tolerance = TOL)
})

test_that("test trace_canon", {
  X <- Constant(rbind(c(1.0, 5.0), c(9.0, 14.0)))
  Y <- matrix_trace(X)
  canon <- Dgp2Dcp.trace_canon(Y, Y@args)
  canon_matrix <- canon[[1]]
  constraints <- canon[[2]]
  expect_equal(length(constraints), 0)
  expect_true(is_scalar(canon_matrix))
  expected <- log(exp(1.0), exp(14.0))
  expect_equal(expected, value(canon_matrix), tolerance = TOL)
})

test_that("test one_minus_pos", {
  x <- Variable(pos = TRUE)
  obj <- Maximize(x)
  constr <- list(one_minus_pos(x) >= 0.4)
  problem <- Problem(obj, constr)
  result <- solve(problem, gp = TRUE)
  expect_equal(result$value, 0.6, tolerance = TOL)
  expect_equal(result$getValue(x), 0.6, tolerance = TOL)
})

test_that("test qp solver not allowed", {
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(x))
  expect_error(solve(problem, solver = "OSQP", gp = TRUE))
})

# TODO: Enable test once sum_largest is implemented.
# test_that("test paper example sum_largest", {
#   x <- Variable(4, pos = TRUE)
#   x1 <- x[1]
#   x2 <- x[2]
#   x3 <- x[3]
#   x4 <- x[4]
#   obj <- Minimize(sum_largest(hstack(list(3*x1^(0.5) * x2^(0.5), x1*x2 + 0.5*x2*x4^3, x3)), 2))
#   constr <- list(x1*x2*x3 >= 16)
#   p <- Problem(obj, constr)
#   solve(p, gp = TRUE)   # Smoke test.
# })

test_that("test paper example one_minus_pos", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x*y)
  constr <- list((y * one_minus_pos(x/y))^2 >= 1, x >= y/3)
  problem <- Problem(obj, constr)
  solve(problem, gp = TRUE)   # Smoke test.
})

test_that("test paper example eye_minus_inv", {
  X <- Variable(2,2, pos = TRUE)
  obj <- Minimize(matrix_trace(eye_minus_inv(X)))
  constr <- list(geo_mean(diag(X)) == 0.1)
  problem <- Problem(obj, constr)
  solve(problem, gp = TRUE, solver = "SCS")   # Smoke test.
})

test_that("test paper example exp_log", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x*y)
  constr <- list(exp(y/x) <= log(y))
  problem <- Problem(obj, constr)
  solve(problem, gp = TRUE)   # Smoke test.
})

test_that("test pf matrix completion", {
  X <- Variable(3,3, pos = TRUE)
  obj <- Minimize(pf_eigenvalue(X))
  known_indices <- cbind(c(1,1,2,3,3), c(1,3,2,1,2))
  constr <- list(X[known_indices] == c(1.0, 1.9, 0.8, 3.2, 5.9), 
                 X[1,2] * X[2,1] * X[2,3] * X[3,3] == 1.0)
  problem <- Problem(obj, constr)
  solve(problem, gp = TRUE)   # Smoke test.
})

test_that("test rank one nmf", {
  X <- Variable(3,3, pos = TRUE)
  x <- Variable(3, pos = TRUE)
  y <- Variable(3, pos = TRUE)
  xy <- vstack(list(x[1]*y, x[2]*y, x[3]*y))
  R <- max_elemwise(multiply(X, (xy)^(-1.0)), multiply(X^(-1.0), xy))
  objective <- sum(R)
  constraints <- list(X[1,1] == 1.0, X[1,3] == 1.9, X[2,2] == 0.8, X[3,1] == 3.2, X[3,2] == 5.9,
                      x[1] * x[2] * x[3] == 1.0)
  # Smoke test.
  prob <- Problem(Minimize(objective), constraints)
  solve(prob, gp = TRUE)
})

test_that("test documentation problem", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  
  objective_fn <- x*y*z
  constraints <- list(4*x*y*z + 2*x*z <= 10, x <= 2*y, y <= 2*x, z >= 1)
  problem <- Problem(Maximize(objective_fn), constraints)
  solve(problem, gp = TRUE)   # Smoke test.
})

test_that("test solver error", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod <- x*y
  dgp <- Problem(Minimize(prod))
  dgp2dcp <- Dgp2Dcp()
  inverse_data <- perform(dgp2dcp, dgp)[[2]]
  soln <- Solution(SOLVER_ERROR, NA, list(), list(), list())
  dgp_soln <- invert(dgp2dcp, soln, inverse_data)
  expect_equal(dgp_soln$status, SOLVER_ERROR)
})
