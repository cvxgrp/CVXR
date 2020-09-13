context("test-g01-lin_ops")

NEG <- "neg"
PARAM <- "param"
DENSE_CONST <- "dense_const"
SPARSE_CONST <- "sparse_const"
SCALAR_CONST <- "scalar_const"
SUM_ENTRIES <- "sum_entries"
VARIABLE <- "variable"

create_const <- CVXR:::create_const
create_eq <- CVXR:::create_eq
create_leq <- CVXR:::create_leq
create_param <- CVXR:::create_param
create_var <- CVXR:::create_var
get_expr_vars <- CVXR:::get_expr_vars

lo.neg_expr <- CVXR:::lo.neg_expr
lo.sum_expr <- CVXR:::lo.sum_expr
lo.sum_entries <- CVXR:::lo.sum_entries

test_that("Test creating a variable", {
  skip_on_cran()
  var <- create_var(c(5, 4), var_id = 1)
  expect_equal(var$dim, c(5, 4))
  expect_equal(var$data, 1)
  expect_equal(length(var$args), 0)
  expect_equal(var$type, VARIABLE)
})

test_that("Test creating a parameter", {
  skip_on_cran()
  A <- Parameter(5, 4)
  var <- create_param(A, c(5, 4))
  expect_equal(var$dim, c(5, 4))
  expect_equal(length(var$args), 0)
  expect_equal(var$type, PARAM)
})

test_that("Test creating a constant", {
  skip_on_cran()
  # Scalar constant
  dim <- c(1, 1)
  mat <- create_const(1.0, dim)
  expect_equal(mat$dim, dim)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, SCALAR_CONST)
  expect_equal(mat$data, 1.0)

  # Dense matrix constant
  dim <- c(5, 4)
  mat <- create_const(matrix(1, nrow = dim[1], ncol = dim[2]), dim)
  expect_equal(mat$dim, dim)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, DENSE_CONST)
  expect_equal(mat$data, matrix(1, nrow = dim[1], ncol = dim[2]))

  # Sparse matrix constant
  dim <- c(5, 5)
  mat <- create_const(Matrix::sparseMatrix(i = 1:5, j = 1:5, x = 1), dim, sparse = TRUE)
  expect_equal(mat$dim, dim)
  expect_equal(length(mat$args), 0)
  expect_equal(mat$type, SPARSE_CONST)
  expect_equivalent(as.matrix(mat$data), diag(rep(1, 5)))
})

test_that("Test adding lin expr", {
  skip_on_cran()
  dim <- c(5, 4)
  x <- create_var(dim)
  y <- create_var(dim)

  # Expanding dict.
  add_expr <- lo.sum_expr(list(x, y))
  expect_equal(add_expr$dim, dim)
  expect_equal(length(add_expr$args), 2)
})

test_that("Test getting vars from an expression", {
  skip_on_cran()
  dim <- c(5, 4)
  x <- create_var(dim)
  y <- create_var(dim)
  A <- create_const(matrix(1, nrow = dim[1], ncol = dim[2]), dim)

  # Expanding dict.
  add_expr <- lo.sum_expr(list(x, y, A))
  vars_ <- get_expr_vars(add_expr)
  expect_equal(vars_[[1]], list(x$data, dim))
  expect_equal(vars_[[2]], list(y$data, dim))
})

test_that("Test negating an expression", {
  skip_on_cran()
  dim <- c(5, 4)
  var <- create_var(dim)
  expr <- lo.neg_expr(var)
  expect_equal(length(expr$args), 1)
  expect_equal(expr$dim, dim)
  expect_equal(expr$type, NEG)
})

test_that("Test creating an equality constraint", {
  skip_on_cran()
  dim <- c(5, 5)
  x <- create_var(dim)
  y <- create_var(dim)
  lh_expr <- lo.sum_expr(list(x, y))
  value <- matrix(1, nrow = dim[1], ncol = dim[2])
  rh_expr <- create_const(value, dim)
  constr <- create_eq(lh_expr, rh_expr)
  expect_equal(constr$dim, dim)
  vars_ <- get_expr_vars(constr$expr)
  ref <- list(list(x$data, dim), list(y$data, dim))
  expect_equal(vars_, ref)
})

test_that("Test creating a less than or equal constraint", {
  skip_on_cran()
  dim <- c(5, 5)
  x <- create_var(dim)
  y <- create_var(dim)
  lh_expr <- lo.sum_expr(list(x, y))
  value <- matrix(1, nrow = dim[1], ncol = dim[2])
  rh_expr <- create_const(value, dim)
  constr <- create_leq(lh_expr, rh_expr)
  expect_equal(constr$dim, dim)
  vars_ <- get_expr_vars(constr$expr)
  ref <- list(list(x$data, dim), list(y$data, dim))
  expect_equal(vars_, ref)
})

test_that("Test sum entries op", {
  skip_on_cran()
  dim <- c(5, 5)
  x <- create_var(dim)
  expr <- lo.sum_entries(x, c(1, 1))
  expect_equal(expr$dim, c(1, 1))
  expect_equal(length(expr$args), 1)
  expect_equal(expr$type, SUM_ENTRIES)
})
