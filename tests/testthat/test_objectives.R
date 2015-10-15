x <- Variable(name = "x")
y <- Variable(3, name = "y")
z <- Variable(name = "z")

test_that("test the Minimize class", {
  exp <- x + z
  obj <- Minimize(exp)
  
  # canon <- canonicalize(obj)
  # new_obj <- canon[[1]]
  # constraints <- canon[[2]]
  
  # expect_equal(length(constraints), 0)
  # expect_error(canonicalize(Minimize(y)))
})

test_that("test the Maximize class", {
  exp <- x + z
  obj <- Maximize(exp)
  
  # canon <- canonicalize(obj)
  # new_obj <- canon[[1]]
  # constraints <- canon[[2]]
  
  # expect_equal(length(constraints), 0)
  # expert_error(canonicalize(Maximize(y)))
})

test_that("test is_dcp for Minimize and Maximize", {
  expect_true(is_dcp(Minimize(NormInf(x))))
  expect_false(is_dcp(Minimize(-NormInf(x))))
  
  expect_false(is_dcp(Maximize(NormInf(x))))
  expect_true(is_dcp(Maximize(-NormInf(x))))
})