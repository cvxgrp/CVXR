context("test-g01-semidefinite_vars")
TOL <- 1e-6

X <- Variable(2, 2, PSD = TRUE)
Y <- Variable(2, 2)
Fmat <- rbind(c(1,0), c(0,-1))

test_that("test that results are symmetric", {
  skip_on_cran()
  M <- Variable(3, 3, PSD = TRUE)
  C1 <- rbind(c(0, 0, 1/2), c(0, 0, 0), c(1/2, 0, 1))
  C2 <- rbind(c(0, 0, 0), c(0, 0, 1/2), c(0, 1/2, 1))
  x1 <- Variable(3, 3, PSD = TRUE)
  x2 <- Variable(3, 3, PSD = TRUE)
  constraints <- list(M + C1 == x1)
  constraints <- c(constraints, M + C2 == x2)
  objective <- Minimize(matrix_trace(M))
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  M_val <- result$getValue(M)
  expect_equal(M_val, t(M_val), tolerance = TOL)
})

test_that("SDP in objective and constraint", {
  skip_on_cran()
  # PSD in objective.
  obj <- Minimize(sum((X - Fmat)^2))
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = 1e-4)

  Xres <- result$getValue(X)
  expect_equal(Xres[1,1], 1, tolerance = 1e-3)
  expect_equal(Xres[1,2], 0, tolerance = TOL)
  expect_equal(Xres[2,1], 0, tolerance = TOL)
  expect_equal(Xres[2,2], 0, tolerance = TOL)

  # PSD in constraint.
  # ECHU: note to self, apparently this is a source of redundancy.
  obj <- Minimize(sum((Y - Fmat)^2))
  p <- Problem(obj, list(Y == Variable(2, 2, PSD = TRUE)))
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = 1e-2)

  Yres <- result$getValue(Y)
  expect_equal(Yres[1,1], 1, tolerance = 1e-3)
  expect_equal(Yres[1,2], 0, tolerance = TOL)
  expect_equal(Yres[2,1], 0, tolerance = TOL)
  expect_equal(Yres[2,2], 0, tolerance = 1e-3)

  # Index into semidef
  obj <- Minimize((X[1,1] - 1)^2 + (X[2,1] - 2)^2 + (X[2,2] - 4)^2)
  p <- Problem(obj, list())
  result <- solve(p)
  print(result$getValue(X))
  expect_equal(result$value, 0, tolerance = 1e-5)

  Xres <- result$getValue(X)
  expect_equal(Xres[1,1], 1, tolerance = 1e-2)
  expect_equal(Xres[1,2], 2, tolerance = 1e-2)
  expect_equal(Xres[2,1], 2, tolerance = 1e-2)
  expect_equal(Xres[2,2], 4, tolerance = 1e-3)
})
