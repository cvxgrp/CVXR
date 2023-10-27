context("test-g01-interfaces")

test_that("Test sign for a given interface", {
  mat <- rbind(c(1,2,3,4), c(3,4,5,6))
  expect_equal(intf_sign(mat), c(TRUE, FALSE))    # Positive.
  expect_equal(intf_sign(-mat), c(FALSE, TRUE))   # Negative.
  expect_equal(intf_sign(0*mat), c(TRUE, TRUE))   # Zero.

  mat <- rbind(c(-1,2,3,4), c(3,4,5,6))
  expect_equal(intf_sign(mat), c(FALSE, FALSE))   # Unknown.
})

test_that("Test R vector interface", {
  # Const to matrix.
  mat <- c(1,2,3)
  expect_equal(intf_dim(mat), c(3,1))
  mat <- c(1,2)
  expect_equal(intf_dim(mat), c(2,1))
  mat <- matrix(2, nrow = 4, ncol = 3)
  expect_equal(intf_dim(mat), c(4,3))
  expect_equal(intf_index(mat, c(2,3)), 2)
  
  # Reshape.
  mat <- rbind(c(1,2,3), c(3,4,5))
  mat <- matrix(mat, nrow = 6, ncol = 1, byrow = TRUE)
  expect_equal(intf_index(mat, c(5,1)), 4)
  
  # Index.
  mat <- rbind(c(1,2,3,4), c(3,4,5,6))
  expect_equal(intf_index(mat, c(1,2)), 3)
  mat <- intf_index(mat, list(seq(1,4,2), seq(1,2)))
  expect_equal(as.vector(mat), c(2,4,4,6))
  
  # Scalars and matrices.
  scalar <- 2
  mat <- c(1,2,3)
  expect_true(all(scalar*mat == c(2,4,6)))
  expect_true(all(scalar - mat == c(1,0,-1)))
  
  # Dimensions.
  expect_equal(intf_dim(c(1,2,3)), c(3,1))
})

test_that("test R matrix interface", {
  # Const to matrix.
  mat <- matrix(c(1,2,3))
  expect_equal(dim(mat), c(3,1))
  mat <- matrix(c(1,2,3), nrow = 3)
  expect_equal(mat[1,1], 1)
  mat <- intf_scalar_matrix(2, c(4,3))
  expect_equal(intf_dim(mat), c(4,3))
  expect_equal(intf_index(mat, c(1,2)), 2)
  
  # Reshape.
  mat <- rbind(c(1,2,3), c(3,4,5))
  mat <- intf_reshape(mat, c(6,1))
  expect_equal(intf_index(mat, c(5,1)), 4)
  
  # Index.
  mat <- rbind(c(1,2,3,4), c(3,4,5,6))
  expect_equal(intf_index(mat, c(1,2)), 3)
  mat <- intf_index(mat, list(seq(1,4,2), seq(1,2)))
  expect_false(any(mat - rbind(c(2,4), c(4,6))))
})

test_that("Test Matrix sparse interface", {
  if(!require(Matrix)) {
    print("Matrix library not found. Skipping test.")
    return()
  }
  
  # Const to matrix.
  mat <- Matrix(c(1,2,3), sparse = TRUE)
  expect_equal(intf_dim(mat), c(3,1))
  # C <- sparseMatrix(i = c(1, 2, 3, 1, 1), j = c(1, 1, 1, 2, 3), x = c(1, 1, 1, 1, 1))
  # expect_equal(intf_dim(mat), c(3,3))
  
  # Identity.
  cmp_mat <- diag(4)
  mat <- Matrix(cmp_mat, sparse = TRUE)
  expect_equal(intf_dim(mat), intf_dim(cmp_mat))
  expect_equal(sum((mat - cmp_mat) != 0), 0)
  
  # Scalar matrix.
  mat <- Matrix(matrix(2, nrow = 4, ncol = 3), sparse = TRUE)
  expect_equal(intf_dim(mat), c(4,3))
  expect_equal(intf_index(mat, c(2,3)), 2)
  
  # Reshape.
  mat <- rbind(c(1,2,3), c(3,4,5))
  mat <- intf_reshape(mat, c(6,1))
  expect_equal(intf_index(mat, c(5,1)), 4)
  
  # Test scalars.
  scalar <- intf_scalar_matrix(1, c(1,3))
  expect_equal(dim(scalar), c(1,3))
  
  # Index.
  mat <- Matrix(rbind(c(1,2,3,4), c(3,4,5,6)), sparse = TRUE)
  expect_equal(intf_index(mat, c(1,2)), 3)
  mat <- intf_index(mat, list(seq(1,4,2), seq(1,2)))
  expect_false(any(mat - rbind(c(2,4), c(4,6))))
  
  # Scalar value.
  mat <- Matrix(diag(1), sparse = TRUE)
  expect_equal(intf_scalar_value(mat), 1.0)
})
