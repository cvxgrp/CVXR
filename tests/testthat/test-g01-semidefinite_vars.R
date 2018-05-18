TOL <- 1e-6

X <- Semidef(2)
Y <- Variable(2, 2)
Fmat <- rbind(c(1,0), c(0,-1))

test_that("SDP in objective and constraint", {
  # SDP in objective
  obj <- Minimize(sum((X - Fmat)^2))
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = 1e-4)
  
  Xres <- result$getValue(X)
  expect_equal(Xres[1,1], 1, tolerance = 1e-3)
  expect_equal(Xres[1,2], 0, tolerance = TOL)
  expect_equal(Xres[2,1], 0, tolerance = TOL)
  expect_equal(Xres[2,2], 0, tolerance = TOL)
  
  # SDP in constraint
  obj <- Minimize(sum((Y - Fmat)^2))
  p <- Problem(obj, list(Y == Semidef(2)))
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
