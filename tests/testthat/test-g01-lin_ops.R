##require(Matrix)

test_that("Test creating a variable", {
  var <- create_var(c(5, 4), var_id = 1)
  expect_equal(var$size, c(5, 4))
  expect_equal(var$data, 1)
  expect_equal(length(var$args), 0)
  expect_equal(var$type, VARIABLE)
})

test_that("Test creating a parameter", {
  A <- Parameter(5, 4)
  var <- create_param(A, c(5, 4))
  expect_equal(var$size, c(5, 4))
  expect_equal(length(var$args), 0)
  expect_equal(var$type, PARAM)
})

test_that("Test creating a constant", {
  # Scalar constant
  size <- c(1, 1)
  mat <- create_const(1.0, size)
  expect_equal(mat$size, size)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, SCALAR_CONST)
  expect_equal(mat$data, 1.0)

  # Dense matrix constant
  size <- c(5, 4)
  mat <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)
  expect_equal(mat$size, size)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, DENSE_CONST)
  expect_equal(mat$data, matrix(1, nrow = size[1], ncol = size[2]))

  # Sparse matrix constant
  size <- c(5, 5)
  mat <- create_const(Matrix::sparseMatrix(i = 1:5, j = 1:5, x = 1), size, sparse = TRUE)
  expect_equal(mat$size, size)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, SPARSE_CONST)
  expect_equivalent(as.matrix(mat$data), diag(rep(1, 5)))
})

test_that("Test adding lin expr", {
  size <- c(5, 4)
  x <- create_var(size)
  y <- create_var(size)

  # Expanding dict.
  add_expr <- lo.sum_expr(list(x, y))
  expect_equal(add_expr$size, size)
  expect_equal(length(add_expr$args), 2)
})

test_that("Test getting vars from an expression", {
  size <- c(5, 4)
  x <- create_var(size)
  y <- create_var(size)
  A <- create_const(matrix(1, nrow = size[1], ncol = size[2]), size)

  # Expanding dict.
  add_expr <- lo.sum_expr(list(x, y, A))
  vars_ <- get_expr_vars(add_expr)
  expect_equal(vars_[[1]], list(x$data, size))
  expect_equal(vars_[[2]], list(y$data, size))
})

test_that("Test negating an expression", {
  size <- c(5, 4)
  var <- create_var(size)
  expr <- lo.neg_expr(var)
  expect_equal(length(expr$args), 1)
  expect_equal(expr$size, size)
  expect_equal(expr$type, NEG)
})

test_that("Test creating an equality constraint", {
  size <- c(5, 5)
  x <- create_var(size)
  y <- create_var(size)
  lh_expr <- lo.sum_expr(list(x, y))
  value <- matrix(1, nrow = size[1], ncol = size[2])
  rh_expr <- create_const(value, size)
  constr <- create_eq(lh_expr, rh_expr)
  expect_equal(constr$size, size)
  vars_ <- get_expr_vars(constr$expr)
  ref <- list(list(x$data, size), list(y$data, size))
  expect_equal(vars_, ref)
})

test_that("Test creating a less than or equal constraint", {
  size <- c(5, 5)
  x <- create_var(size)
  y <- create_var(size)
  lh_expr <- lo.sum_expr(list(x, y))
  value <- matrix(1, nrow = size[1], ncol = size[2])
  rh_expr <- create_const(value, size)
  constr <- create_leq(lh_expr, rh_expr)
  expect_equal(constr$size, size)
  vars_ <- get_expr_vars(constr$expr)
  ref <- list(list(x$data, size), list(y$data, size))
  expect_equal(vars_, ref)
})

test_that("Test sum entries op", {
  size <- c(5, 5)
  x <- create_var(size)
  expr <- lo.sum_entries(x)
  expect_equal(expr$size, c(1, 1))
  expect_equal(length(expr$args), 1)
  expect_equal(expr$type, SUM_ENTRIES)
})
