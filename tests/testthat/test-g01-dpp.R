context("test-g01-dpp")
TOL <- 1e-6
SOLVER <- "ECOS"

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
test_that("test basic equality constraint", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  dgp <- Problem(Minimize(x), list(x == alpha))
  
  expect_true(is_dgp(dgp@objective, dpp = TRUE))
  expect_true(is_dgp(dgp@constraints[[1]], dpp = TRUE))
  expect_true(is_dgp(dgp, dpp = TRUE))
  dgp2dcp <- Dgp2Dcp(dgp)
  
  dcp <- reduce(dgp2dcp)
  expect_true(is_dpp(dcp))
  
  result <- solve(dgp, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  
  value(alpha) <- 2.0
  dgp <- Problem(Minimize(x), list(x == alpha))
  result <- solve(dgp, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(x), 2.0, tolerance = TOL)
})

test_that("test basic inequality constraint", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  constraint <- list(x + alpha <= x)
  expect_true(is_dgp(consraint[[1]], dpp = TRUE))
  expect_true(is_dgp(Problem(Minimize(1), constraint), dpp = TRUE))
})

test_that("test nonlla equality constraint not dpp", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  constraint <- list(x == x + alpha)
  expect_false(is_dgp(constraint[[1]], dpp = TRUE))
  expect_false(is_dgp(Problem(Minimize(1), constraint), dpp = TRUE))
})

test_that("test nonllcvx inequality constraint not dpp", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  constraint <- list(x <= x + alpha)
  expect_false(is_dgp(constraint[[1]], dpp = TRUE))
  expect_false(is_dgp(Problem(Minimize(1), constraint), dpp = TRUE))
})

test_that("test param monomial is dpp", {
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)
  
  monomial <- alpha^1.2 * beta^0.5 * kappa^3 * kappa^2
  expect_true(is_dgp(monomial, dpp = TRUE))
})

test_that("test param posynomial is dpp", {
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)
  
  monomial <- alpha^1.2 * beta^0.5 * kappa^3 * kappa^2
  posynomial <- monomial + alpha^2 * beta^3
  expect_true(is_dgp(posynomial, dpp = TRUE))
})

test_that("test mixed monomial is dpp", {
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)
  tau <- Variable(pos = TRUE)
  
  monomial <- alpha^1.2 * beta^0.5 * kappa^3 * kappa^2 * tau
  expect_true(is_dgp(monomial, dpp = TRUE))
})

test_that("test mixed posynomial is dpp", {
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)
  tau <- Variable(pos = TRUE)
  
  monomial <- alpha^1.2 * beta^0.5 * kappa^3 * kappa^2 * tau
  posynomial <- (monomial + monomial)^3
  expect_true(is_dgp(posynomial, dpp = TRUE))
})

test_that("test nested power not dpp", {
  alpha <- Parameter(value = 1.0)
  x <- Variable(pos = TRUE)
  
  pow1 <- x^alpha
  expect_true(is_dgp(pow1, dpp = TRUE))
  
  pow2 <- pow1^alpha
  expect_false(is_dgp(pow2, dpp = TRUE))
})

test_that("test non dpp problem raises error", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  dgp <- Problem(Minimize((alpha*x)^alpha), list(x == alpha))
  expect_true(is_dgp(dgp@objective))
  expect_false(is_dgp(dgp@objective, dpp = TRUE))
  
  expect_error(solve(dgp, solver = "SCS", gp = TRUE, enforce_dpp = TRUE))
  
  result <- solve(dgp, solver = "SCS", gp = TRUE, enforce_dpp = FALSE)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
})

test_that("test basic monomial", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  beta <- Parameter(pos = TRUE, value = 2.0)
  x <- Variable(pos = TRUE)
  monomial <- alpha*beta*x
  problem <- Problem(Minimize(monomial), list(x == alpha))
  
  expect_true(is_dgp(problem))
  expect_true(is_dgp(problem, dpp = TRUE))
  expect_false(is_dpp(problem, "dcp"))
  
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$value, 2.0, tolerance = TOL)
  
  value(alpha) <- 3.0
  # problem <- Problem(Minimize(monomial), list(x == alpha))
  expect_equal(result$getValue(x), 3.0, tolerance = TOL)
  # 3*2*3 == 18
  expect_equal(result$value, 18.0, tolerance = TOL)
})

test_that("test basic posynomial", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  beta <- Parameter(pos = TRUE, value = 2.0)
  kappa <- Parameter(pos = TRUE, value = 3.0)
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  monomial_one <- alpha*beta*x
  monomial_two <- beta*kappa*x*y
  posynomial <- monomial_one + monomial_two
  problem <- Problem(Minimize(posynomial), list(x == alpha, y == beta))
  
  expect_true(is_dgp(problem))
  expect_true(is_dgp(problem, dpp = TRUE))
  expect_false(is_dpp(problem, "dcp"))
  
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 2.0, tolerance = TOL)
  # 1*2*1 + 2*3*1*2 == 2 + 12 == 14.
  expect_equal(result$value, 14.0, tolerance = 1e-3)
  
  value(alpha) <- 4.0
  value(beta) <- 5.0
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(x), 4.0, tolerance = TOL)
  expect_equal(result$getValue(y), 5.0, tolerance = TOL)
  # 4*5*4 + 5*3*4*5 == 80 + 300 == 380.
  expect_equal(result$value, 380.0, tolerance = 1e-3)
})

test_that("test basic gp", {
  xyz <- Variable(3, pos = TRUE)
  x <- xyz[1]
  y <- xyz[2]
  z <- xyz[3]
  a <- Parameter(pos = TRUE, value = 2.0)
  b <- Parameter(pos = TRUE, value = 1.0)
  constraints <- list(a*x*y + a*x*z + a*y*z <= b, x >= a*y)
  problem <- Problem(Minimize(1/(x*y*z)), constraints)
  expect_true(is_dgp(problem, dpp = TRUE))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 15.59, tolerance = 1e-2)
})

test_that("test maximum", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  alpha <- Parameter(value = 0.5)
  beta <- Parameter(pos = TRUE, value = 3.0)
  kappa <- Parameter(pos = TRUE, value = 1.0)
  tau <- Parameter(pos = TRUE, value = 4.0)
  
  prod1 <- x*y^alpha
  prod2 <- beta*x*y^alpha
  obj <- Minimize(max_elemwise(prod1, prod2))
  constr <- list(x == kappa, y == tau)
  
  problem <- Problem(obj, constr)
  expect_true(is_dgp(problem, dpp = TRUE))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # max(1*2, 3*1*2) = 6.
  expect_equal(result$value, 6.0, tolerance = 1e-4)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
  
  value(alpha) <- 2
  value(beta) <- 0.5
  value(kappa) <- 2.0   # x
  value(tau) <- 3.0     # y
  # obj <- Minimize(max_elemwise(prod1, prod2))
  # constr <- list(x == kappa, y == tau)
  # problem <- Problem(obj, constr)
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # max(2*9, 0.5*2*9) == 18.
  expect_equal(result$value, 18.0, tolerance = 1e-4)
  expect_equal(result$getValue(x), 2.0, tolerance = TOL)
  expect_equal(result$getValue(y), 3.0, tolerance = TOL)
})

test_that("test max", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  alpha <- Parameter(value = 0.5)
  beta <- Parameter(pos = TRUE, value = 3.0)
  kappa <- Parameter(pos = TRUE, value = 1.0)
  tau <- Parameter(pos = TRUE, value = 4.0)
  
  prod1 <- x*y^alpha
  prod2 <- beta*x*y^alpha
  obj <- Minimize(max(hstack(prod1, prod2)))
  constr <- list(x == kappa, y == tau)
  
  problem <- Problem(obj, constr)
  expect_true(is_dgp(problem, dpp = TRUE))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # max(1*2, 3*1*2) = 6.
  expect_equal(result$value, 6.0, tolerance = 1e-4)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
  
  value(alpha) <- 2
  value(beta) <- 0.5
  value(kappa) <- 2.0   # x
  value(tau) <- 3.0     # y
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # max(2*9, 0.5*2*9) == 18.
  expect_equal(result$value, 18.0, tolerance = 1e-4)
  expect_equal(result$getValue(x), 2.0, tolerance = TOL)
  expect_equal(result$getValue(y), 3.0, tolerance = TOL)
})

test_that("test param in exponent and elsewhere", {
  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(x^alpha), list(x == alpha))
  
  expect_true(is_dgp(problem, dpp = TRUE))
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 1.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  
  # Re-solve (which goes through a separate code path).
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 27.0, tolerance = TOL)
  expect_equal(result$getValue(x), 3.0, tolerance = TOL)
})

test_that("test minimum", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  beta <- Parameter(pos = TRUE, value = 3.0, name = "beta")
  prod1 <- x*y^alpha
  prod2 <- beta*x*y^alpha
  posy <- prod1 + prod2
  obj <- Maximize(min_elemwise(prod1, prod2, 1/posy))
  constr <- list(x == alpha, y == 4.0)
  
  dgp <- Problem(obj, constr)
  result <- solve(dgp, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # prod1 = 1*4, prod2 = 3*4 = 12, 1/posy = 1/(3 + 12).
  expect_equal(result$value, 1.0/(4.0 + 12.0), tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
  
  value(alpha) <- 2.0
  # prod1 <- x*y^alpha
  # prod2 <- beta*x*y^alpha
  # posy <- prod1 + prod2
  # obj <- Maximize(min_elemwise(prod1, prod2, 1/posy))
  # constr <- list(x == alpha, y == 4.0)
  # dgp <- Problem(obj, constr)
  # prod1 = 2*16, prod2 = 3*2*16 = 96, 1/posy = 1/(32 + 96).
  result <- solve(dgp, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 1.0/(32.0 + 96.0), tolerance = TOL)
  expect_equal(result$getValue(x), 2.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
})

test_that("test min", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  beta <- Parameter(pos = TRUE, value = 3.0, name = "beta")
  prod1 <- x*y^alpha
  prod2 <- beta*x*y^alpha
  posy <- prod1 + prod2
  obj <- Maximize(min_entries(hstack(prod1, prod2, 1/posy)))
  constr <- list(x == alpha, y == 4.0)
  
  dgp <- Problem(obj, constr)
  result <- solve(dgp, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  # prod1 = 1*4, prod2 = 3*4 = 12, 1/posy = 1/(3 + 12).
  expect_equal(result$value, 1.0/(4.0 + 12.0), tolerance = TOL)
  expect_equal(result$getValue(x), 1.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
  
  value(alpha) <- 2.0
  # prod1 <- x*y^alpha
  # prod2 <- beta*x*y^alpha
  # posy <- prod1 + prod2
  # obj <- Maximize(min_entries(hstack(prod1, prod2, 1/posy)))
  # constr <- list(x == alpha, y == 4.0)
  # dgp <- Problem(obj, constr)
  # prod1 = 2*16, prod2 = 3*2*16 = 96, 1/posy = 1/(32 + 96).
  result <- solve(dgp, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 1.0/(32.0 + 96.0), tolerance = TOL)
  expect_equal(result$getValue(x), 2.0, tolerance = TOL)
  expect_equal(result$getValue(y), 4.0, tolerance = TOL)
})

test_that("test div", {
  alpha <- Parameter(pos = TRUE, value = 3.0)
  beta <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  
  p <- Problem(Minimize(x*y), list(y/alpha <= x, y >= beta))
  result <- solve(p, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 1.0/3.0, tolerance = TOL)
  expect_equal(result$getValue(x), 1.0/3.0, tolerance = TOL)
  expect_equal(result$getValue(y), 1.0, tolerance = TOL)
  
  value(beta) <- 2.0
  p <- Problem(Minimize(x*y), list(y/alpha <= x, y >= beta))
  result <- solve(p, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 4.0/3.0, tolerance = TOL)
  expect_equal(result$getValue(x), 2.0/3.0, tolerance = TOL)
  expect_equal(result$getValue(y), 2.0, tolerance = TOL)
})

test_that("test one minus pos", {
  x <- Variable(pos = TRUE)
  obj <- Maximize(x)
  alpha <- Parameter(pos = TRUE, value = 0.1)
  constr <- list(one_minus_pos(alpha + x) >= 0.4)
  problem <- Problem(obj, constr)
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 0.5, tolerance = TOL)
  expect_equal(result$getValue(x), 0.5, tolerance = TOL)
  
  value(alpha) <- 0.4
  # constr <- list(one_minus_pos(alpha + x) >= 0.4)
  # problem <- Problem(obj, constr)
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 0.2, tolerance = TOL)
  expect_equal(result$getValue(x), 0.2, tolerance = TOL)
})

test_that("test pf matrix completion", {
  X <- Variable(3, 3, pos = TRUE)
  obj <- Minimize(pf_eigenvalue(X))
  known_indices <- cbind(c(1,1,2,3,3), c(1,3,2,1,2))
  constr <- list(X[known_indices] == c(1.0, 1.9, 0.8, 3.2, 5.9),
                 X[1,2] * X[2,1] * X[2,3] * X[3,3] == 1.0)
  problem <- Problem(obj, constr)
  # Smoke test.
  result <- solve(problem, SOLVER, gp = TRUE)
  optimal_value <- result$value
  
  param <- Parameter(length(known_values), pos = TRUE, value = 0.5*known_values)
  constr <- list(X[known_indices] == param,
                 X[1,2] * X[2,1] * X[2,3] * X[3,3] == 1.0)
  problem <- Problem(obj, constr)
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  
  # Now change param to point to known_value, and check we recover the correct optimal value.
  value(param) <- known_values
  # constr <- list(X[known_indices] == param,
  #                X[1,2] * X[2,1] * X[2,3] * X[3,3] == 1.0)
  # problem <- Problem(obj, constr)
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, optimal_value, tolerance = TOL)
})

test_that("test rank one nmf", {
  X <- Variable(3, 3, pos = TRUE)
  x <- Variable(3, pos = TRUE)
  y <- Variable(3, pos = TRUE)
  xy <- vstack(x[1]*y, x[2]*y, x[3]*y)
  R <- max_elemwise(mul_elemwise(X, (xy)^(-1.0)), 
                    mul_elemwise(X^(-1.0), xy))
  objective <- sum(R)
  constraints <- list(X[1,1] == 1.0, X[1,3] == 1.9, X[2,2] == 0.8, 
                      X[3,1] == 3.2, X[3,2] == 5.9, x[1]*x[2]*x[3] == 1.0)
  # Smoke test.
  prob <- Problem(Minimize(objective), constraints)
  result <- solve(prob, SOLVER, gp = TRUE)
  optimal_value <- result$value
  
  param <- Parameter(value = -2.0)
  R <- max_elemwise(mul_elemwise(X, (xy)^param),
                    mul_elemwise(X^param, xy))
  objective <- sum(R)
  prob <- Problem(Minimize(objective), constraints)
  result <- solve(prob, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  
  # Now change param to point to known_value, and check we recover the correct
  # optimal value.
  value(param) <- -1.0
  result <- solve(prob, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, optimal_value, tolerance = TOL)
})

test_that("test documentation prob", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 4.0)
  b <- Parameter(pos = TRUE, value = 2.0)
  c <- Parameter(pos = TRUE, value = 10.0)
  d <- Parameter(pos = TRUE, value = 1.0)
  
  objective_fn <- x*y*z
  constraints <- list(a*x*y*z + b*x*z <= c, x <= b*y, y <= b*x, z >= d)
  problem <- Problem(Maximize(objective_fn), constraints)
  # Smoke test.
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
})

test_that("test sum scalar", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  w <- Variable(pos = TRUE)
  h <- Variable(pos = TRUE)
  problem <- Problem(Minimize(h), list(w*h >= 8, sum(alpha + w) <= 5))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(h), 2, tolerance = TOL)
  expect_equal(result$getValue(w), 4, tolerance = TOL)
  
  value(alpha) <- 4.0
  problem <- Problem(Minimize(h), list(w*h >= 8, sum(alpha + w) <= 5))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 8, tolerance = TOL)
  expect_equal(result$getValue(h), 8, tolerance = TOL)
  expect_equal(result$getValue(w), 1, tolerance = TOL)
})

test_that("test sum vector", {
  alpha <- Parameter(2, pos = TRUE, value = c(1.0, 1.0))
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  problem <- Problem(Minimize(sum(h)), list(mul_elemwise(w, h) >= 20, sum(alpha + w) <= 10))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(h), matrix(c(5, 5)), tolerance = 1e-3)
  expect_equal(result$getValue(w), matrix(c(4, 4)), tolerance = 1e-3)
  
  value(alpha) <- c(4.0, 4.0)
  problem <- Problem(Minimize(sum(h)), list(mul_elemwise(w, h) >= 20, sum(alpha + w) <= 10))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 40, tolerance = TOL)
  expect_equal(result$getValue(h), matrix(c(20, 20)), tolerance = 1e-3)
  expect_equal(result$getValue(w), matrix(c(1, 1)), tolerance = 1e-3)
})

test_that("test sum squares vector", {
  alpha <- Parameter(2, pos = TRUE, value = c(1.0, 1.0))
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  problem <- Problem(Minimize(sum_squares(alpha + h)), list(mul_elemwise(w, h) >= 20, sum(alpha + w) <= 10))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(w), matrix(c(4, 4)), tolerance = 1e-3)
  expect_equal(result$getValue(h), matrix(c(5, 5)), tolerance = 1e-3)
  expect_equal(result$value, 6^2 + 6^2, tolerance = TOL)
  
  value(alpha) <- c(4.0, 4.0)
  problem <- Problem(Minimize(sum_squares(alpha + h)), list(mul_elemwise(w, h) >= 20, sum(alpha + w) <= 10))
  result <- solve(problem, SOLVEr, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$getValue(w), matrix(c(1, 1)), tolerance = 1e-3)
  expect_equal(result$getValue(h), matrix(c(20, 20)), tolerance = 1e-3)
  expect_equal(result$value, 24^2 + 24^2, tolerance = 1e-3)
})

test_that("test sum matrix", {
  w <- Variable(2, 2, pos = TRUE)
  h <- Variable(2, 2, pos = TRUE)
  alpha <- Parameter(pos = TRUE, value = 1.0)
  problem <- Problem(Minimize(alpha*sum(h)), list(mul_elemwise(w, h) >= 10, sum(w) <= 20))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 8, tolerance = 1e-4)
  expect_equal(result$getValue(h), rbind(c(2, 2), c(2, 2)), tolerance = 1e-4)
  expect_equal(result$getValue(w), rbind(c(5, 5), c(5, 5)), tolerance = 1e-4)
  
  value(alpha) <- 2.0
  result <- solve(problem, SOLVER, gp = TRUe, enforce_dpp = TRUE)
  expect_equal(result$value, 16, tolerance = 1e-4)
  
  w <- Variable(2, 2, pos = TRUE)
  h <- Parameter(2, 2, pos = TRUE)
  value(h) <- matrix(1, nrow = 2, ncol = 2)
  value(alpha) <- 1.0
  problem <- Problem(Minimize(sum(alpha*h)), list(w == h))
  result <- solve(problem, SOLVEr, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 4.0, tolerance = 1e-4)
  
  value(h) <- 2.0*matrix(1, nrow = 2, ncol = 2)
  # problem <- Problem(Minimize(sum(alpha*h)), list(w == h))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 8.0, tolerance = 1e-4)
  
  value(h) <- 3.0*matrix(1, nrow = 2, ncol = 2)
  # problem <- Problem(Minimize(sum(alpha*h)), list(w == h))
  result <- solve(problem, SOLVER, gp = TRUE, enforce_dpp = TRUE)
  expect_equal(result$value, 12.0, tolerance = 1e-4)
})

test_that("test exp", {
  x <- Variable(4, pos = TRUE)
  c <- Parameter(4, pos = TRUE)
  expr <- exp(mul_elemwise(c, x))
  expect_true(is_dgp(expr, dpp = TRUE))
  
  expr <- exp(t(c) %*% x)
  expect_true(is_dgp(expr, dpp = TRUE))
})

test_that("test log", {
  x <- Variable(4, pos = TRUE)
  c <- Parameter(4, pos = TRUE)
  expr <- log(mul_elemwise(c, x))
  expect_true(is_dgp(expr, dpp = TRUE))
  
  expr <- log(t(c) %*% x)
  expect_false(is_dgp(expr, dpp = TRUE))
})

test_that("test gmatmul", {
  x <- Variable(2, pos = TRUE)
  A <- Parameter(2, 2)
  value(A) <- rbind(c(-5, 2), c(1, -3))
  b <- matrix(c(3, 2))
  expr <- gmatmul(A, x)
  problem <- Problem(Minimize(1.0), list(expr == b))
  expect_true(is_dgp(problem, dpp = TRUE))
  result <- solve(problem, solver = "SCS", gp = TRUE, enforce_dpp = TRUE)
  sltn <- exp(base::solve(value(A), log(b)))
  expect_equal(result$getValue(x), sltn, tolerance = TOL)
  
  x_par <- Parameter(2, pos = TRUE)
  expr <- gmatmul(A, x_par)
  expect_false(is_dgp(expr, dpp = TRUE))
  expect_true(is_dgp(expr, dpp = FALSE))
})
