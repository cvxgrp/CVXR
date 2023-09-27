context("test-g01-suppfunc")
TOL <- 1e-6

# Test the implementation of support function atoms.
# Relevant source code includes:
#    SuppFuncAtom in atoms.R
#    SuppFunc in transforms.R
#    Dcp2Cone.suppfunc_canon in dcp2cone

test_that("test Rn", {
  set.seed(1)
  n <- 5
  x <- Variable(n)
  sigma <- suppfunc(x, list())
  a <- matrix(rnorm(n), nrow = n, ncol = 1)
  y <- Variable(n)
  cons <- list(sigma(y - a) <= 0)   # "<= num" for any num >= 0 is valid.
  objective <- Minimize(t(a) %*% y)
  prob <- Problem(objective, cons)
  result <- solve(prob, solver = "ECOS")
  actual <- result$value
  expected <- t(a) %*% a
  expect_lte(abs(actual - expected), 1e-6)
  actual <- value(y)
  expected <- a
  expect_lte(base::norm(actual - expected, "2"), 1e-6)
  viol <- result$getViolation(cons[[1]])
  expect_lte(viol, 1e-8)
})

test_that("test vector1norm", {
  n <- 3
  set.seed(1)
  a <- matrix(rnorm(n), nrow = n, ncol = 1)
  x <- Variable(n)
  sigma <- suppfunc(x, list(norm(x - a, 1) <= 1))
  y <- matrix(rnorm(n), nrow = n, ncol = 1)
  y_var <- Variable(n)
  prob <- Problem(Minimize(sigma(y_var)), list(y == y_var))
  result <- solve(prob, solver = "ECOS")
  actual <- result$value
  expected <- t(a) %*% y + base::norm(y, "I")
  expect_lte(abs(actual - expected), 1e-5)
  expect_lte(abs(result$getValue(prob@objective@expr) - result$value), 1e-5)
})

test_that("test vector2norm", {
  n <- 3
  set.seed(1)
  a <- matrix(rnorm(n), nrow = n, ncol = 1)
  x <- Variable(n)
  sigma <- suppfunc(x, list(norm(x - a, 2) <= 1))
  y <- matrix(rnorm(n), nrow = n, ncol = 1)
  y_var <- Variable(n)
  prob <- Problem(Minimize(sigma(y_var)), list(y == y_var))
  result <- solve(prob, solver = "ECOS")
  actual <- result$value
  expected <- t(a) %*% y + base::norm(y, "2")
  expect_lte(abs(actual - expected), 1e-5)
  expect_lte(abs(result$getValue(prob@objective@expr) - result$value), 1e-5)
})

test_that("test rectangular variable", {
  set.seed(2)
  rows <- 4
  cols <- 2
  a <- matrix(rnorm(rows*cols), nrow = rows, ncol = cols)
  x <- Variable(rows, cols)
  sigma <- suppfunc(x, list(x[,1] == 0))
  y <- Variable(rows, cols)
  cons <- list(sigma(y - a) <= 0)
  objective <- Minimize(sum_squares(flatten(y)))
  prob <- Problem(objective, cons)
  result <- solve(prob, solver = "ECOS")
  expect <- cbind(matrix(0, nrow = rows, ncol = 1), a[,2])
  actual <- result$getValue(y)
  expect_lte(base::norm(actual - expect, "2"), 1e-6)
  viol <- result$getViolation(cons[[1]])
  expect_lte(viol, 1e-6)
})

test_that("test psd dualcone", {
  set.seed(5)
  n <- 3
  X <- Variable(n, n)
  sigma <- suppfunc(X, list(X %>>% 0))
  A <- matrix(rnorm(n^2), nrow = n, ncol = n)
  Y <- Variable(n, n)
  objective <- Minimize(norm(as.vector(A) + flatten(Y)))
  cons <- list(sigma(Y) <= 0)   # Y is negative definite.
  prob <- Problem(objective, cons)
  result <- solve(prob, solver = "SCS", eps = 1e-8)
  viol <- result$getViolation(cons[[1]])
  expect_lte(viol, 1e-6)
  eigs <- eigen(result$getValue(Y))$values
  expect_lte(max(eigs), 1e-6)
})

test_that("test largest singvalue", {
  set.seed(3)
  rows <- 3
  cols <- 4
  A <- matrix(rnorm(rows*cols), nrow = rows, ncol = cols)
  A_sv <- svd(A)
  X <- Variable(rows, cols)
  sigma <- suppfunc(X, list(sigma_max(X) <= 1))
  Y <- Variable(rows, cols)
  cons <- list(Y == A)
  prob <- Problem(Minimize(sigma(Y)), cons)
  result <- solve(prob, solver = "SCS", eps = 1e-8)
  actual <- result$value
  expect <- sum(A_sv$d)
  expect_lte(abs(actual - expect), 1e-6)
})

test_that("test expcone 1", {
  x <- Variable(1)
  tempcons <- list(exp(x[1]) <= exp(1), exp(-x[1]) <= exp(1))
  sigma <- suppfunc(x, tempcons)
  y <- Variable(1)
  obj_expr <- y[1]
  cons <- list(sigma(y) <= 1)   # This just means -1 <= y[1] <= 1.
  prob <- Problem(Minimize(obj_expr), cons)
  result <- solve(prob, solver = "ECOS")
  viol <- result$getViolation(cons[[1]])
  expect_lte(viol, 1e-6)
  expect_lte(abs(result$getValue(y) - (-1)), 1e-6)
})

test_that("test expcone 2", {
  x <- Variable(3)
  tempcons <- list(sum(x) <= 1.0, sum(x) >= 0.1, x >= 0.01, kl_div(x[2], x[1]) + x[2] - x[1] + x[3] <= 0)
  sigma <- suppfunc(x, tempcons)
  y <- Variable(3)
  a <- matrix(c(-3, -2, -1))   # This is negative of objective in mosek_conif example.
  expr <- -sigma(y)
  objective <- Maximize(expr)
  cons <- list(y == a)
  prob <- Problem(objective, cons)
  result <- solve(prob, solver = "ECOS")
  # Check fo expected objective value.
  epi_actual <- result$value
  direct_actual <- result$getValue(expr)
  expect <- 0.235348211
  expect_lte(abs(epi_actual - expect), 1e-6)
  expect_lte(abs(direct_actual - expect), 1e-6)
})

test_that("test basic lmi", {
  set.seed(4)
  A <- matrix(rnorm(n^2), nrow = n, ncol = n)
  A <- t(A) %*% A
  X <- Variable(n, n)   # Will fail if you try PSD = TRUE or symmetric = TRUE.
  sigma <- suppfunc(X, list(0 %<<% X, lambda_max(X) <= 1))
  Y <- Variable(n, n)
  cons <- list(Y == A)
  expr <- sigma(Y)
  prob <- Problem(Minimize(expr), cons)   # opt value of support func would be at X = I.
  result <- solve(prob, solver = "SCS", eps = 1e-8)
  actual1 <- result$value            # Computed with epigraph.
  actual2 <- result$getValue(expr)   # Computed by evaluating support function, as a maximization problem.
  expect_lte(abs(actual1 - actual2), 1e-6)
  expect <- sum(diag(A))
  expect_lte(abs(actual1 - expect), 1e-4)
})

test_that("test invalid solver", {
  n <- 3
  x <- Variable(n)
  sigma <- suppfunc(x, list(norm(x - rnorm(n), 2) <= 1))
  y_var <- Variable(n)
  prob <- Problem(Minimize(sigma(y_var)), list(rnorm(n) == y_var))
  expect_error(solve(prob, solver = "OSQP"), ".*could not be reduced to a QP.*")
})

test_that("test invalid variable", {
  x <- Variable(2, 2, symmetric = TRUE)
  expect_error(suppfunc(x, list()))
})

test_that("test invalid constraint", {
  x <- Variable(3)
  a <- Parameter(3)
  cons <- list(t(a) %*% x == 1)
  expect_error(suppfunc(x, cons))
})
