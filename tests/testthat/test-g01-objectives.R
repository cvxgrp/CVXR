context("test-g01-objectives")

x <- Variable(name = "x")
y <- Variable(3, name = "y")
z <- Variable(name = "z")

test_that("test string representations", {
  obj <- Minimize(x)
  expect_equal(as.character(obj), paste("Minimize(", as.character(x), ")", sep = ""))
  obj <- Minimize(2*x)
  expect_equal(as.character(obj), paste("Minimize(", as.character(2*x), ")", sep = ""))
  
  obj <- Maximize(x)
  expect_equal(as.character(obj), paste("Maximize(", as.character(x), ")", sep = ""))
  obj <- Maximize(2*x)
  expect_equal(as.character(obj), paste("Maximize(", as.character(2*x), ")", sep = ""))
})

test_that("test the Minimize class", {
  expr <- x + z
  obj <- Minimize(expr)
  expect_equal(as.character(obj), paste("minimize", name(expr)))
  
  canon <- canonical_form(obj)
  new_obj <- canon[[1]]
  constraints <- canon[[2]]
  # expect_equal(name(constraints[[1]]), name(new_obj == expr))
  
  # For affine objectives, there should be no constraints.
  expect_equal(length(constraints), 0)
  expect_error(canonical_form(Minimize(y)), 
               "The 'minimize' objective must resolve to a scalar.", fixed = TRUE)
  
  # Test copy with args = NULL.
  copy <- copy(obj)
  expect_true(identical(class(copy), class(obj)))

  # A new object is constructed, so copy@args == obj@args, but copy@args is not obj@args.
  expect_equal(copy@args, obj@args)
  expect_false(identical(copy@args, obj@args))
  
  # Test copy with new args.
  copy <- copy(obj, args = list(square(z)))
  expect_true(identical(class(copy), class(obj)))
  expect_true(identical(copy@args[[1]]@args[[1]], z))
})

test_that("test the Maximize class", {
  expr <- x + z
  obj <- Maximize(expr)
  expect_equal(as.character(obj), paste("maximize", name(expr)))
  
  canon <- canonical_form(obj)
  new_obj <- canon[[1]]
  constraints <- canon[[2]]
  # expect_equal(name(constraints[[1]]), name(new_obj == expr))
  
  # For affine objectives, there should be no constraints.
  expect_equal(length(constraints), 0)
  expect_error(canonical_form(Maximize(y)),
               "The 'maximize' objective must resolve to a scalar.", fixed = TRUE)
  
  # Test copy with args = NULL.
  copy <- copy(obj)
  expect_true(identical(class(copy), class(obj)))
  
  # A new object is constructed, so copy@args == obj@args, but copy@args is not obj@args.
  expect_equal(copy@args, obj@args)
  expect_false(identical(copy@args, obj@args))
  
  # Test copy with new args.
  copy <- copy(obj, args = list(-square(x)))
  expect_true(identical(class(copy), class(obj)))
  expect_true(identical(copy@args[[1]]@args[[1]], x))
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
  
  # Addition.
  expect_true(is_dcp(Minimize(expr1) + Minimize(expr2)))
  expect_true(is_dcp(Maximize(-expr1) + Maximize(-expr2)))
  
  # Test Minimize + Maximize.
  expect_error(Minimize(expr1) + Maximize(-expr2), "Problem does not follow DCP rules.", fixed = TRUE)
  expect_true(is_dcp(Minimize(expr1) - Maximize(-expr2)))
  
  # Multiplication (alpha is a positive scalar).
  expect_true(is_dcp(alpha*Minimize(expr1)))
  expect_true(is_dcp(alpha*Maximize(-expr1)))
  expect_true(is_dcp(-alpha*Maximize(-expr1)))
  expect_true(is_dcp(-alpha*Maximize(-expr1)))
})
