context("test-g01-dpp")
TOL <- 1e-6

test_that("test multiply scalar params not dpp", {
  x <- Parameter()
  product <- x*x
  expect_false(is_dpp(product))
  expect_true(is_dcp(product))
})

test_that("test matmul params not dpp", {
  X <- Parameter(4, 4)
  product <- X %*% X
  expect_true(is_dcp(product))
  expect_false(is_dpp(product))
})

test_that("test multiply param and variable is dpp", {
  x <- Parameter()
  y <- Variable()
  product <- x*y
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

test_that("test multiply variable and param is dpp", {
  x <- Parameter()
  y <- Variable()
  product <- mul_elemwise(y, x)
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

test_that("test multiply nonlinear param and variable is not dpp", {
  x <- Parameter()
  y <- Variable()
  product <- exp(x)*y
  expect_false(is_dpp(product))
})

test_that("test multiply nonlinear nonneg param and nonneg variable is not dpp", {
  x <- Parameter(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  product <- exp(x)*y
  expect_false(is_dpp(product))
  expect_true(is_dcp(product))
})

test_that("test multiply affine param and variable is dpp", {
  x <- Parameter()
  y <- Variable()
  product <- (x + x)*y
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

test_that("test multiply param plus var times const", {
  x <- Parameter()
  y <- Variable()
  product <- (x + y)*5
  expect_true(is_convex(product))
  expect_true(is_dcp(product))
  expect_true(is_dpp(product))
})

test_that("test multiply param and nonlinear variable is dpp", {
  x <- Parameter(nonneg = TRUE)
  y <- Variable()
  product <- x*exp(y)
  expect_true(is_convex(product))
  expect_true(is_dcp(product))
  expect_true(is_dpp(product))
})

test_that("test nonlinear equality not dpp", {
  x <- Variable()
  a <- Parameter()
  constraint <- list(x == norm(a))
  expect_false(is_dcp(constraint[[1]], dpp = TRUE))
  problem <- Problem(Minimize(0), constraint)
  expect_false(is_dcp(problem, dpp = TRUE))
})

test_that("test nonconvex inequality not dpp", {
  x <- Variable()
  a <- Parameter()
  constraint <- list(x <= norm(a))
  expect_false(is_dcp(constraints[[1]], dpp = TRUE))
  problem <- Problem(Minimize(0), constraint)
  expect_false(is_dcp(problem, dpp = TRUE))
})

test_that("test solve multiply param plus var times const", {
  x <- Parameter()
  y <- Variable()
  product <- (x + y)*5
  expect_true(is_dpp(product))
  value(x) <- 2.0
  problem <- Problem(Minimize(product), list(y == 1))
  result <- solve(problem, "SCS")
  expect_equal(result$value, 15, tolerance = TOL)
})

test_that("test paper example is dpp", {
  Fparm <- Parameter(2, 2)
  x <- Variable(2, 1)
  g <- Parameter(2, 1)
  lambd <- Parameter(nonneg = TRUE)
  objective <- norm(Fparm %*% x - g) + lambd*norm(x)
  constraints <- list(x >= 0)
  problem <- Problem(Minimize(objective), constraints)
  expect_true(is_dpp(objective))
  expect_true(is_dpp(constraints[[1]]))
  expect_true(is_dpp(problem))
})

test_that("test non dcp expression is not dpp", {
  x <- Parameter()
  expr <- exp(log(x))
  expect_false(is_dpp(expr))
})

test_that("test can solve non dpp problem", {
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x*x), list(x == y))
  expect_false(is_dpp(problem))
  expect_true(is_dcp(problem))
  expect_equal(solve(problem, "SCS")$value, 25)
  value(x) <- 3
  # problem <- Problem(Minimize(x*x), list(x == y))
  expect_equal(solve(problem, "SCS"), 9)
})

test_that("test chain data for non dpp problem evals params", {
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x*x), list(x == y))
  chain <- get_problem_data(problem, "SCS")[[2]]
  expect_false(is_dpp(problem))
  expect_true("EvalParams" %in% sapply(chain@reductions, class))
})

test_that("test chain data for dpp problem does not eval params", {
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x + y), list(x == y))
  chain <- get_problem_data(problem, "SCS")[[2]]
  expect_false("EvalParams" %in% sapply(chain@reductions, class))
})

test_that("test param quad form not dpp", {
  x <- Variable(2, 1)
  P <- Parameter(2, 2, PSD = TRUE)
  value(P) <- diag(2)
  y <- quad_form(x, P)
  expect_false(is_dpp(y))
  expect_true(is_dcp(y))
})

test_that("test const quad form is dpp", {
  x <- Variable(2, 1)
  P <- diag(2)
  y <- quad_form(x, P)
  expect_true(is_dpp(y))
  expect_true(is_dcp(y))
})

test_that("test paper example logreg is dpp", {
  N <- 3
  n <- 2
  beta <- Variable(n, 1)
  b <- Variable(1, 1)
  X <- Parameter(N, n)
  Y <- matrix(1, nrow = N, ncol = 1)
  lambd1 <- Parameter(nonneg = TRUE)
  lambd2 <- Parameter(nonneg = TRUE)
  log_likelihood <- (1/N)*sum(mul_elemwise(Y, X %*% beta + b) - 
                              t(log_sum_exp(t(hstack(matrix(0, nrow = N, ncol = 1), X %*% beta + b)),
                                          axis = 2, keepdims = TRUE)))
  regularization <- -lambd1*norm(beta, 1) - lambd2*sum_squares(beta)
  problem <- Problem(Maximize(log_likelihood + regularization))
  expect_true(is_dpp(log_likelihood))
  expect_true(is_dcp(problem))
  expect_true(is_dpp(problem))
})

test_that("test paper example stoch control", {
  n <- 3
  m <- 3
  x <- Parameter(n, 1)
  P_sqrt <- Parameter(m, m)
  P_21 <- Parameter(n, m)
  q <- Parameter(m, 1)
  u <- Variable(m, 1)
  y <- Variable(n, 1)
  objective <- 0.5*sum_squares(P_sqrt %*% u) + t(x) %*% y + t(q) %*% u
  problem <- Problem(Minimize(objective), list(norm(u) <= 0.5, y == P_21 %*% u))
  expect_true(is_dpp(problem))
  expect_true(is_dcp(problem))
})

test_that("test paper example relu", {
  n <- 2
  x <- Parameter(n)
  y <- Variable(n)
  objective <- Minimize(sum_squares(y - x))
  constraints <- list(y >= 0)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))
  value(x) <- c(5, 5)
  result <- solve(problem, "SCS", eps = 1e-8)
  expect_equal(result$getValue(y), result$getValue(x), tolerance = TOL)
  value(x) <- c(-4, -4)
  # problem <- Problem(objective, constraints)
  result <- solve(problem, "SCS", eps = 1e-8)
  expect_equal(result$getValue(y), matrix(c(0, 0)), tolerance = TOL)
})

test_that("test paper example opt net qp", {
  m <- 3
  n <- 2
  G <- Parameter(m, n)
  h <- Parameter(m, 1)
  p <- Parameter(n, 1)
  y <- Variable(n, 1)
  objective <- Minimize(0.5*sum_squares(y - p))
  constraints <- list(G %*% y <= h)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))
})

test_that("test paper example ellipsoidal constraints", {
  n <- 2
  A_sqrt <- Parameter(n, n)
  z <- Parameter(n)
  p <- Parameter(n)
  y <- Variable(n)
  slack <- new("Variable", dim = dim(y))
  objective <- Minimize(0.5*sum_squares(y - p))
  constraints <- list(0.5*sum_squares(A_sqrt %*% slack) <= 1, slack == y - z)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))
})

test_that("test non dpp powers", {
  s <- Parameter(1, nonneg = TRUE)
  x <- Variable(1)
  obj <- Maximize(x + s)
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 2, tolerance = 1e-3)
  
  s <- Parameter(1, nonneg = TRUE)
  x <- Variable(1)
  obj <- Maximize(x + s^2)
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 2, tolerance = 1e-3)
  
  s <- Parameter(1, nonneg = TRUE)
  x <- Variable(1)
  obj <- Maximize(mul_elemwise(x, s^2))
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 1, tolerance = 1e-3)
})

test_that("test ignore dpp", {
  # Test the ignore_dpp flag.
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x + y), list(x == y))
  expect_true(is_dpp(problem))
  expect_true(is_dcp(problem))
  # Basic solve functionality.
  result <- solve(problem, "SCS", ignore_dpp = TRUE)
  expect_equal(result$value, 10, tolerance = TOL)
  
  # enforce_dpp clashes with ignore_dpp.
  expect_error(solve(problem, "SCS", enforce_dpp = TRUE, ignore_dpp = TRUE))
})

###########################
#                         #
# DGP Tests with DPP Flag #
#                         #
###########################


