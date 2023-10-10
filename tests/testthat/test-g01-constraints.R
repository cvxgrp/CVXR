context("test-g01-constraints")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

SOC <- CVXR:::SOC
save_value <- CVXR:::save_value

test_that("test the EqConstraint class", {
  skip_on_cran()
  constr <- x == z
  expect_equal(name(constr), "x == z")
  expect_equal(dim(constr), c(2,1))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))

  x <- save_value(x, 2)
  z <- save_value(z, 2)
  constr <- x == z
  expect_true(constr_value(constr))
  x <- save_value(x, 3)
  constr <- x == z
  expect_false(constr_value(constr))

  value(x) <- c(2,1)
  value(z) <- c(2,2)
  constr <- x == z
  expect_false(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,1)), tolerance = TOL)

  value(z) <- c(2,1)
  constr <- x == z
  expect_true(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,0)))
  expect_equal(residual(constr), matrix(c(0,0)))

  expect_error(x == y)
})

test_that("test the IneqConstraint class", {
  skip_on_cran()
  constr <- x <= z
  expect_equal(name(constr), "x <= z")
  expect_equal(dim(constr), c(2,1))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  x <- save_value(x, 1)
  z <- save_value(z, 2)
  constr <- x <= z
  expect_true(constr_value(constr))
  x <- save_value(x, 3)
  constr <- x <= z
  expect_false(constr_value(constr))

  value(x) <- c(2,1)
  value(z) <- c(2,0)
  constr <- x <= z
  expect_false(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,1)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,1)), tolerance = TOL)

  value(z) <- c(2,2)
  constr <- x <= z
  expect_true(constr_value(constr))
  expect_equal(violation(constr), matrix(c(0,0)), tolerance = TOL)
  expect_equal(residual(constr), matrix(c(0,0)), tolerance = TOL)
  
  # Incompatible dimensions.
  expect_error(x <= y)
})

test_that("Test the PSD constraint %>>%", {
  skip_on_cran()
  constr <- A %>>% B
  expect_equal(name(constr), "A + -B >> 0")
  expect_equal(dim(constr), c(2,2))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  A <- save_value(A, rbind(c(2,-1), c(1,2)))
  B <- save_value(B, rbind(c(1,0), c(0,1)))
  constr <- A %>>% B
  expect_true(constr_value(constr))
  expect_equal(violation(constr), 0, tolerance = TOL)
  expect_equal(residual(constr), 0, tolerance = TOL)

  B <- save_value(B, rbind(c(3,0), c(0,3)))
  constr <- A %>>% B
  expect_false(constr_value(constr))
  expect_equal(violation(constr), 1, tolerance = TOL)
  expect_equal(residual(constr), 1, tolerance = TOL)

  expect_error(x %>>% 0, "Non-square matrix in positive definite constraint.")
})

test_that("Test the PSD constraint %<<%", {
  skip_on_cran()
  constr <- A %<<% B
  expect_equal(name(constr), "B + -A >> 0")
  expect_equal(dim(constr), c(2,2))

  # Test value and dual_value
  expect_true(is.na(dual_value(constr)))
  expect_error(constr_value(constr))
  B <- save_value(B, rbind(c(2,-1), c(1,2)))
  A <- save_value(A, rbind(c(1,0), c(0,1)))
  constr <- A %<<% B
  expect_true(constr_value(constr))
  A <- save_value(A, rbind(c(3,0), c(0,3)))
  constr <- A %<<% B
  expect_false(constr_value(constr))

  expect_error(x %<<% 0, "Non-square matrix in positive definite constraint.")
})

test_that("test the >= operator", {
  skip_on_cran()
  constr <- z >= x
  expect_equal(name(constr), "x <= z")
  expect_equal(dim(constr), c(2,1))
  
  # Incompatible dimensions.
  expect_error(y >= x)
})

test_that("test the SOC class", {
  skip_on_cran()
  exp <- x + z
  scalar_exp <- a + b
  constr <- SOC(scalar_exp, exp)
  expect_equal(size(constr), 3)
  
  # Test invalid dimensions.
  expect_error(SOC(Variable(1), Variable(1, 4), axis = 2), 
               "Argument dimensions (1,1) and (1,4), with axis = 2, are incompatible", fixed = TRUE)
  
  # Test residual.
  # 1D.
  n <- 5
  x0 <- 0:(n-1)
  t0 <- 2
  xv <- Variable(n, value = x0)
  tv <- Variable(value = t0)
  resid <- residual(SOC(tv, xv))
  expect_equal(ndim(resid), 0)
  dist <- sum_squares(xv - x0) + square(tv - t0)
  prob <- Problem(Minimize(dist), list(SOC(tv, xv)))
  result <- solve(prob)
  expect_equal(sqrt(result$getValue(dist)), resid, tolerance = TOL)
  
  # 2D, axis = 1.
  n <- 5
  k <- 3
  x0 <- matrix(0:(n*k-1), nrow = n, ncol = k, byrow = TRUE)
  t0 <- matrix(c(1,2,3))
  xv <- Variable(n, k, value = x0)
  tv <- Variable(k, value = t0)
  resid <- residual(SOC(tv, xv, axis = 1))
  expect_equal(dim(resid), c(k,1))
  for(i in 1:k) {
    dist <- sum_squares(xv[i,] - x0[i,]) + sum_squares(tv[i] - t0[i])
    prob <- Problem(Minimize(dist), list(SOC(tv[i], xv[i,])))
    result <- solve(prob)
    expect_equal(sqrt(result$getValue(dist)), resid[i], tolerance = TOL)
  }
  
  # 2D, axis = 2.
  n <- 5
  k <- 3
  x0 <- matrix(0:(n*k-1), nrow = n, ncol = k, byrow = TRUE)
  t0 <- matrix(c(1,2,3))
  xv <- Variable(n, k, value = x0)
  tv <- Variable(k, value = t0)
  resid <- residual(SOC(tv, xv, axis = 2))
  expect_equal(dim(resid), c(k,1))
  for(i in 1:k) {
    dist <- sum_squares(xv[,i] - x0[,i]) + sum_squares(tv[i] - t0[i])
    prob <- Problem(Minimize(dist), list(SOC(tv[i], xv[,i])))
    result <- solve(prob)
    expect_equal(sqrt(result$getValue(dist)), resid[i], tolerance = TOL)
  }
  
  # Test all three cases:
  # 1. t >= ||x||
  # 2. -||x|| < t < ||x||
  # 3. t <= -||x||
  k <- 3
  n <- 3
  x0 <- matrix(1, nrow = k, ncol = n)
  norms <- norm(x0, "2")
  t0 <- matrix(c(2, 0.5, -2))*norms
  xv <- Variable(k, n, value = x0)
  tv <- Variable(k, value = t0)
  resid <- residual(SOC(tv, xv, axis = 1))
  expect_equal(dim(resid), c(k,1))
  for(i in 1:k) {
    dist <- sum_squares(xv[i,] - x0[i,]) + sum_squares(tv[i] - t0[i])
    prob <- Problem(Minimize(dist), list(SOC(t[i], x[i,])))
    result <- solve(prob)
    expect_equal(sqrt(result$getValue(dist)), resid[i], tolerance = 1e-4)
  }
})

test_that("test Pow3D constraint", {
  n <- 3
  set.seed(0)
  alpha <- 0.275
  xv <- Variable(n)
  yv <- Variable(n)
  zv <- Variable(n)
  con <- PowCone3D(xv, yv, zv, alpha)
  
  # Check violation against feasible values.
  x0 <- 0.1 + runif(n)
  y0 <- 0.1 + runif(n)
  z0 <- x0^alpha + y0^(1 - alpha)
  z0[2] <- z0[2]*-1
  
  value(xv) <- x0
  value(yv) <- y0
  value(zv) <- z0
  # con <- PowCone3D(xv, yv, zv, alpha)
  viol <- residual(con)
  expect_gte(viol, 0.99*abs(x1[1]))
  
  # Check invalid constraint data.
  expect_error(con <- PowCone3D(xv, yv, zv, 1.001))
  expect_error(con <- PowCone3D(xv, yv, zv, -0.00001))
})

test_that("test PowND constraint", {
  n <- 4
  Wv <- Variable(n)
  zv <- Variable()
  set.seed(0)
  alpha <- 0.5 + runif(n)
  alpha <- alpha/sum(alpha)
  
  expect_error(con <- PowConeND(Wv, zv, alpha + 0.01))   # Entries don't sum to one.
  expect_error(con <- PowConeND(Wv, zv, matrix(alpha, nrow = n, ncol = 1, byrow = TRUE)))
  expect_error(con <- PowConeND(reshape_expr(Wv, c(n,1)), zv, matrix(alpha, nrow = n, ncol = 1, byrow = TRUE), axis = 1))

  # Compute a violation.
  con <- PowConeND(Wv, zv, alpha)
  W0 <- 0.01 + runif(n)
  z0 <- prod(W0^alpha) + 0.05
  value(W) <- W0
  value(z) <- z0
  # con <- PowConeND(Wv, zv, alpha)
  viol <- violation(con)
  expect_gte(viol, 0.01)
  expect_lte(viol, 0.06)
})

test_that("test chained constraints", {
  # Tests that chaining constraints raises an error.
  error_str <- "Cannot evaluate the truth value of a constraint or chain constraints, e.g., 1 >= x >= 0"
  expect_error(z <= x <= 1, error_str, fixed = TRUE)
  expect_error(x == z == 1, error_str, fixed = TRUE)
  expect_error(as.logical(z <= x), error_str, fixed = TRUE)
})

test_that("test the NonPosConstraint for correctness", {
  n <- 3
  xv <- Variable(n)
  cc <- 0:(n-1)
  prob <- Problem(Maximize(sum(xv)), list(NonPosConstraint(xv - cc)))
  
  # Solve through cone program path.
  result1 <- solve(prob, solver = "ECOS")
  expect_equal(result1$getValue(xv), cc, tolerance = TOL)
  
  # Solve through QP path.
  result2 <- solve(prob, solver = "OSQP")
  expect_equal(result2$getValue(xv), cc, tolerance = TOL)
})

test_that("test dual variables work for NonPosConstraint", {
  n <- 3
  xv <- Variable(n)
  cc <- 0:(n-1)
  prob <- Problem(Maximize(sum(xv)), list((xv - cc) <= 0))
  result <- solve(prob, solver = "ECOS")
  dual <- result$getDualValue(prob@constraints[[1]])
  
  prob <- Problem(Maximize(sum(xv)), list(NonPosConstraint(xv - cc)))
  
  # Solve through cone program path.
  result1 <- solve(prob, solver = "ECOS")
  expect_equal(result1$getDualValue(prob@constraints[[1]]), dual, tolerance = TOL)
  
  # Solve through QP path.
  result2 <- solve(prob, solver = "OSQP")
  expect_equal(result2$getDualValue(prob@constraints[[1]]), dual, tolerance = TOL)
})
