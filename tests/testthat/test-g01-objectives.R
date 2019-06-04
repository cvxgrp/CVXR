context("test-g01-objectives")

x <- Variable(name = "x")
y <- Variable(3, name = "y")
z <- Variable(name = "z")

test_that("test the Minimize class", {
  expr <- x + z
  obj <- Minimize(expr)
  
  canon <- canonicalize(obj)
  new_obj <- canon[[1]]
  constraints <- canon[[2]]
  
  expect_equal(length(constraints), 0)
  expect_error(canonical_form(Minimize(y)))
})

test_that("test the Maximize class", {
  expr <- x + z
  obj <- Maximize(expr)
  
  canon <- canonicalize(obj)
  new_obj <- canon[[1]]
  constraints <- canon[[2]]
  
  expect_equal(length(constraints), 0)
  expect_error(canonical_form(Maximize(y)))
})

test_that("test is_dcp for Minimize and Maximize", {
  expect_true(is_dcp(Minimize(norm_inf(x))))
  expect_false(is_dcp(Minimize(-norm_inf(x))))
  
  expect_false(is_dcp(Maximize(norm_inf(x))))
  expect_true(is_dcp(Maximize(-norm_inf(x))))
})

test_that("test adding objectives", {
  expr1 <- x^2
  expr2 <- x^-1
  alpha <- 2
  
  # Addition
  expect_true(is_dcp(Minimize(expr1) + Minimize(expr2)))
  expect_true(is_dcp(Maximize(-expr1) + Maximize(-expr2)))
  
  # Test Minimize + Maximize
  expect_error(Minimize(expr1) + Maximize(-expr2))
  expect_true(is_dcp(Minimize(expr1) - Maximize(-expr2)))
  
  # Multiplication (alpha is a positive scalar)
  expect_true(is_dcp(alpha*Minimize(expr1)))
  expect_true(is_dcp(alpha*Maximize(-expr1)))
  expect_true(is_dcp(-alpha*Maximize(-expr1)))
  expect_true(is_dcp(-alpha*Maximize(-expr1)))
})
