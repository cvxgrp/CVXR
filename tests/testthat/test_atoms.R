TOL <- 1e-5

a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the NormInf class", {
  exp <- x + y
  atom <- NormInf(exp)
  
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_true(is_convex(atom))
  expect_true(is_concave(-atom))
  expect_equal(curvature(NormInf(atom)), Curvature.CONVEX)
  expect_equal(curvature(NormInf(-atom)), Curvature.CONVEX)
})

test_that("test the Norm1 class", {
  exp <- x + y
  atom <- Norm1(exp)
  
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(curvature(Norm1(atom)), Curvature.CONVEX)
  expect_equal(curvature(Norm1(-atom)), Curvature.CONVEX)
})

test_that("test the Norm2 class", {
  exp <- x + y
  atom <- Norm2(exp)
  
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(curvature(Norm2(atom)), Curvature.CONVEX)
  expect_equal(curvature(Norm2(-atom)), Curvature.CONVEX)
  
  # Test with axis arg
  expr <- Pnorm(A, 2, axis = 1)
  expect_equal(size(expr), c(2, 1))
})

test_that("test the Power class", {
  for(size in list(c(1, 1), c(3, 1), c(2, 3))) {
    x_pow <- Variable(size[1], size[2])
    y_pow <- Variable(size[1], size[2])
    exp <- x_pow + y_pow
    
    for(p in c(0, 1, 2, 3, 2.7, 0.67, -1, -2.3, 4/5)) {
      atom <- Power(exp, p)
      expect_equal(size(atom), size)
      
      if(p > 1 || p < 0)
        expect_equal(curvature(atom), Curvature.CONVEX)
      else if(p == 1)
        expect_equal(curvature(atom), Curvature.AFFINE)
      else if(p == 0)
        expect_equal(curvature(atom), Curvature.CONSTANT)
      else
        expect_equal(curvature(atom), Curvature.CONCAVE)
      
      if(p != 1)
        expect_equal(sign(atom), Sign.POSITIVE)
    }
  }
})

test_that("test the GeoMean class", {
  # TODO: Need to implement fracify for this to work
  # atom <- GeoMean(x)
  # expect_equal(size(atom), c(1, 1))
  # expect_equal(curvature(atom), Curvature.CONCAVE)
  # expect_equal(sign(atom), Sign.POSITIVE)
})

test_that("test the HarmonicMean class", {
  atom <- HarmonicMean(x)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
})

test_that("test the Pnorm class", {
  atom <- Pnorm(x, p = 1.5)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 1)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 2)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = Inf)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 0.5)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = 0.7)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -0.1)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -1)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
  
  atom <- Pnorm(x, p = -1.3)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONCAVE)
  expect_equal(sign(atom), Sign.POSITIVE)
})

test_that("test the QuadOverLin class", {
  atom <- QuadOverLin(Square(x), a)
  expect_equal(curvature(atom), Curvature.CONVEX)
  
  atom <- QuadOverLin(-Square(x), a)
  expect_equal(curvature(atom), Curvature.CONVEX)
  
  atom <- QuadOverLin(Sqrt(x), a)
  expect_equal(curvature(atom), Curvature.UNKNOWN)
  expect_false(is_dcp(atom))
  
  expect_error(QuadOverLin(x, x))
})

test_that("test the arg count for MaxElemwise and MinElemwise", {
  expect_error(MaxElemwise(1))
  expect_error(MinElemwise(1))
})

test_that("test the MatrixFrac class", {
  atom <- MatrixFrac(x, A)
  expect_equal(size(atom), c(1, 1))
  expect_equal(curvature(atom), Curvature.CONVEX)
  
  # Test MatrixFrac size validation
  expect_error(MatrixFrac(x, C))
  expect_error(MatrixFrac(Variable(3), A))
})

test_that("test the sign for MaxEntries", {
  expect_equal(sign(MaxEntries(1)), Sign.POSITIVE)
  expect_equal(sign(MaxEntries(-2)), Sign.NEGATIVE)
  expect_equal(sign(MaxEntries(Variable())), Sign.UNKNOWN)
  expect_equal(sign(MaxEntries(0)), Sign.ZERO)
  
  # Test with axis argument
  expect_equal(size(MaxEntries(Variable(2), axis = 1)), c(2, 1))
  expect_equal(size(MaxEntries(Variable(2), axis = 2)), c(1, 1))
  expect_equal(size(MaxEntries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(size(MaxEntries(Variable(2, 3), axis = 2)), c(1, 3))
  
  # Invalid axis
  expect_error(MaxEntries(x, axis = 4))
})

test_that("test the sign for MinEntries", {
  expect_equal(sign(MinEntries(1)), Sign.POSITIVE)
  expect_equal(sign(MinEntries(-2)), Sign.NEGATIVE)
  expect_equal(sign(MinEntries(Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinEntries(0)), Sign.ZERO)
  
  # Test with axis argument
  expect_equal(size(MinEntries(Variable(2), axis = 1)), c(2, 1))
  expect_equal(size(MinEntries(Variable(2), axis = 2)), c(1, 1))
  expect_equal(size(MinEntries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(size(MinEntries(Variable(2, 3), axis = 2)), c(1, 3))
  
  # Invalid axis
  expect_error(MinEntries(x, axis = 4))
})

test_that("test sign logic for MaxElemwise", {
  expect_equal(sign(MaxElemwise(1, 2)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, Variable())), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, -2)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(1, 0)), Sign.POSITIVE)
  
  expect_equal(sign(MaxElemwise(Variable(), 0)), Sign.POSITIVE)
  expect_equal(sign(MaxElemwise(Variable(), Variable())), Sign.UNKNOWN)
  expect_equal(sign(MaxElemwise(Variable(), -2)), Sign.UNKNOWN)
  
  expect_equal(sign(MaxElemwise(0, 0)), Sign.ZERO)
  expect_equal(sign(MaxElemwise(0, -2)), Sign.ZERO)
  
  expect_equal(sign(MaxElemwise(-3, -2)), Sign.NEGATIVE)
  
  # Many args
  expect_equal(sign(MaxElemwise(-2, Variable(), 0, -1, Variable(), -1)), Sign.POSITIVE)
  
  # Promotion
  expect_equal(sign(MaxElemwise(1, Variable(2))), Sign.POSITIVE)
  expect_equal(size(MaxElemwise(1, Variable(2))), c(2,1))
})

test_that("test sign logic for MinElemwise", {
  expect_equal(sign(MinElemwise(1, 2)), Sign.POSITIVE)
  expect_equal(sign(MinElemwise(1, Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinElemwise(1, -2)), Sign.NEGATIVE)
  expect_equal(sign(MinElemwise(1, 0)), Sign.ZERO)
  
  expect_equal(sign(MinElemwise(Variable(), 0)), Sign.NEGATIVE)
  expect_equal(sign(MinElemwise(Variable(), Variable())), Sign.UNKNOWN)
  expect_equal(sign(MinElemwise(Variable(), -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(0, 0)), Sign.ZERO)
  expect_equal(sign(MinElemwise(0, -2)), Sign.NEGATIVE)
  
  expect_equal(sign(MinElemwise(-3, -2)), Sign.NEGATIVE)
  
  # Many args
  expect_equal(sign(MinElemwise(-2, Variable(), 0, -1, Variable(), 1)), Sign.NEGATIVE)
  
  # Promotion
  expect_equal(sign(MinElemwise(-1, Variable(2))), Sign.NEGATIVE)
  expect_equal(size(MinElemwise(-1, Variable(2))), c(2, 1))
})

test_that("test the SumEntries class", {
  expect_equal(sign(SumEntries(1)), Sign.POSITIVE)
  expect_equal(sign(SumEntries(c(1, -1))), Sign.UNKNOWN)
  expect_equal(curvature(SumEntries(c(1, -1))), Curvature.CONSTANT)
  expect_equal(sign(SumEntries(Variable(2))), Sign.UNKNOWN)
  expect_equal(size(SumEntries(Variable(2))), c(1,1))
  expect_equal(curvature(SumEntries(Variable(2))), Curvature.AFFINE)
  
  # Mixed curvature
  expect_equal(curvature(SumEntries( c(1,-1) %*% Square(Variable(2)) )), Curvature.UNKNOWN)
  
  # Test with axis argument
  expect_equal(size(SumEntries(Variable(2), axis = 1)), c(2, 1))
  expect_equal(size(SumEntries(Variable(2), axis = 2)), c(1, 1))
  expect_equal(size(SumEntries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(size(SumEntries(Variable(2, 3), axis = 2)), c(1, 3))
  
  # Invalid axis
  expect_error(SumEntries(x, axis = 4))
})

test_that("test the MulElemwise class", {
  expect_equal(sign(MulElemwise(c(1, -1), x)), Sign.UNKNOWN)
  expect_equal(curvature(MulElemwise(c(1, -1), x)), Curvature.AFFINE)
  expect_equal(size(MulElemwise(c(1, -1), x)), c(2, 1))
  pos_param <- Parameter(2, sign = "positive")
  neg_param <- Parameter(2, sign = "negative")
  expect_equal(sign(MulElemwise(pos_param, pos_param)), Sign.POSITIVE)
  expect_equal(sign(MulElemwise(pos_param, neg_param)), Sign.NEGATIVE)
  expect_equal(sign(MulElemwise(neg_param, neg_param)), Sign.POSITIVE)
  
  expect_equal(curvature(MulElemwise(neg_param, Square(x))), Curvature.CONCAVE)
  
  # Test promotion
  expect_equal(size(MulElemwise(c(1, -1), 1)), c(2, 1))
  expect_equal(size(MulElemwise(1, C)), size(C))
  expect_error(MulElemwise(x, c(1, -1)))
})

test_that("test the VStack class", {
  atom <- VStack(x, y, x)
  expect_equal(size(atom), c(6, 1))
  
  atom <- VStack(A, C, B)
  expect_equal(size(atom), c(7, 2))
  
  entries <- list()
  for(i in 1:size(x)[1]) {
    for(j in 1:size(x)[2]) {
      entries <- c(entries, x[i, j])
    }
  }
  atom <- VStack(unlist(entries))
  
  expect_error(VStack(C, 1))
})

test_that("test the Reshape class", {
  expr <- Reshape(A, 4, 1)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(4, 1))
  
  expr <- Reshape(expr, 2, 2)
  expect_equal(size(expr), c(2, 2))
  
  expr <- Reshape(Square(x), 1, 2)
  expect_equal(sign(expr), Sign.POSITIVE)
  expect_equal(curvature(expr), Curvature.CONVEX)
  expect_equal(size(expr), c(1, 2))
  
  expect_error(Reshape(C, 5, 4))
})

test_that("test the Vec class", {
  expr <- Vec(C)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(6, 1))
  
  expr <- Vec(x)
  expect_equal(size(expr), c(2, 1))
  
  expr <- Vec(Square(a))
  expect_equal(sign(expr), Sign.POSITIVE)
  expect_equal(curvature(expr), Curvature.CONVEX)
  expect_equal(size(expr), c(1, 1))
})

test_that("test the Diag class", {
  expr <- Diag(x)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(2,2))
  
  expr <- Diag(A)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(2, 1))
  
  expr <- Diag(t(x))
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(2, 2))
  
  expect_error(Diag(C))
})

test_that("test the Trace class", {
  expr <- Trace(A)
  expect_equal(sign(expr), Sign.UNKNOWN)
  expect_equal(curvature(expr), Curvature.AFFINE)
  expect_equal(size(expr), c(1, 1))
  
  expect_error(Trace(C))
})

test_that("test the Log1p class", {
  expr <- Log1p(1)
  expect_equal(sign(expr), Sign.POSITIVE)
  expect_equal(curvature(expr), Curvature.CONSTANT)
  expect_equal(size(expr), c(1, 1))
  expr <- Log1p(-0.5)
  expect_equal(sign(expr), Sign.NEGATIVE)
})

test_that("test the UpperTri class", {
  expect_error(UpperTri(C))
})

test_that("test the Huber class", {
  Huber(x, 1)
  expect_error(Huber(x, -1))
  expect_error(Huber(x, c(1, 1)))
  
  # M parameter
  M <- Parameter(sign = "positive")
  # Valid
  Huber(x, M)
  M@value <- 1
  expect_equal(value(Huber(2, M)), 3, tolerance = TOL)
  # Invalid
  M <- Parameter(sign = "negative")
  expect_error(Huber(x, M))
})

test_that("test the SumLargest class", {
  expect_error(SumLargest(x, -1))
  expect_error(LambdaSumLargest(x, 2.4))
  expect_error(LambdaSumLargest(Variable(2, 2), 2.4))
})

test_that("test the SumSmallest class", {
  expect_error(SumSmallest(x, -1))
  expect_error(LambdaSumSmallest(Variable(2, 2), 2.4))
})

test_that("test the Bmat class", {
  v_np <- matrix(1, nrow = 3, ncol = 1)
  v_00 <- matrix(c(0,0), nrow = 2, ncol = 1)
  v_12 <- matrix(c(1,2), nrow = 2, ncol = 1)
  expr <- Bmat(list(list(v_np, v_np), list(v_00, v_12)))
  expect_equal(size(expr), c(5, 2))
  const <- rbind(cbind(v_np, v_np), cbind(c(0, 0), c(1, 2)))
  expect_equal(value(expr), const)
})

test_that("test the Conv class", {
  a <- matrix(1, nrow = 3, ncol = 1)
  b <- Parameter(2, sign = "positive")
  expr <- Conv(a, b)
  expect_true(is_positive(expr))
  expect_equal(size(expr), c(4, 1))
  b <- Parameter(2, sign = "negative")
  expr <- Conv(a, b)
  expect_true(is_negative(expr))
  expect_error(Conv(x, -1))
  expect_error(Conv(cbind(c(0, 1), c(0, 1)), x))
})

test_that("test the Kron class", {
  a <- matrix(1, nrow = 3, ncol = 2)
  b <- Parameter(2, sign = "positive")
  expr <- Kron(a, b)
  expect_true(is_positive(expr))
  expect_equal(size(expr), c(6, 2))
  b <- Parameter(2, sign = "negative")
  expr <- Kron(a, b)
  expect_true(is_negative(expr))
  expect_error(Kron(x, -1))
})

test_that("test DCP properties of partial optimize", {
  # Evaluate the 1-norm in the usual way (i.e., in epigraph form)
  dims <- 3
  x <- Variable(dims)
  t <- Variable(dims)
  xval <- matrix(rep(-5, dims), nrow = dims, ncol = 1)
  p2 <- Problem(Minimize(SumEntries(t)), list(-t <= x, x <= t))
  # g <- partial_optimize(p2, list(t), list(x))
  # expect_equal(curvature(g), Curvature.CONVEX)
  
  p2 <- Problem(Maximize(SumEntries(t)), list(-t <= x, x <= t))
  # g <- partial_optimize(p2, list(t), list(x))
  # expect_equal(curvature(g), Curvature.CONCAVE)
  
  p2 <- Problem(Maximize(Square(t[1])), list(-t <= x, x <= t))
  # g <- partial_optimize(p2, list(t), list(x))
  # expect_false(is_convex(g))
  # expect_false(is_concave(g))
})

test_that("test the partial_optimize eval 1-norm", {
  # Evaluate the 1-norm in the usual way (i.e., in epigraph form)
  dims <- 3
  x <- Variable(dims)
  t <- Variable(dims)
  xval <- matrix(rep(-5, dims), nrow = dims, ncol = 1)
  p1 <- Problem(Minimize(SumEntries(t)), list(-t <= xval, xval <= t))
  # result1 <- solve(p1)
  
  # Minimize the 1-norm via partial_optimize
  p2 <- Problem(Minimize(SumEntries(t)), list(-t <= x, x <= t))
  # g <- partial_optimize(p2, list(t), list(x))
  # p3 <- Problem(Minimize(g), list(x == xval))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, -result3$optimal_value)
  
  # Try leaving out args
  
  # Minimize the 1-norm via partial_optimize
  p2 <- Problem(Minimize(SumEntries(t)), list(-t <= x, x <= t))
  # g <- partial_optimize(p2, opt_vars = list(t))
  # p3 <- Problem(Minimize(g), list(x == xval))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
  
  # Minimize the 1-norm via partial_optimize
  # g <- partial_optimize(p2, dont_opt_vars = list(x))
  # p3 <- Problem(Minimize(g), list(x == xval))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
  
  # expect_error(partial_optimize(p2))
  # expect_error(partial_optimize(p2, list(), list(x)))
})

test_that("test partial_optimize min 1-norm", {
  # Minimize the 1-norm in the usual way
  dims <- 3
  x <- Variable(dims)
  t <- Variable(dims)
  p1 <- Problem(Minimize(SumEntries(t)), list(-t <= x, x <= t))
  
  # Minimize the 1-norm via partial_optimize
  # g <- partial_optimize(p1, list(t), list(x))
  # p2 <- Problem(Minimize(g))
  # result2 <- solve(p2)
  
  # result1 <- solve(p1)
  # expect_equal(result1$optimal_value, result2$optimal_value)
})

test_that("test partial_optimize simple problem", {
  x <- Variable(1)
  y <- Variable(1)
  
  # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
  p1 <- Problem(Minimize(x+y), list(x+y >= 3, y >= 4, x >= 5))
  # result1 <- solve(p1)
  
  # Solve the two-stage problem via partial_optimize
  p2 <- Problem(Minimize(y), list(x+y >= 3, y >= 4))
  # g <- partial_optimize(p2, list(y), list(x))
  # p3 <- Problem(Minimize(x+g), list(x >= 5))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
})

test_that("test partial_optimize special var", {
  x <- Bool(1)
  y <- Int(1)
  
  # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
  p1 <- Problem(Minimize(x+y), list(x+y >= 3, y >= 4, x >= 5))
  # result1 <- solve(p1)
  
  # Solve the two-stage problem via partial_optimize
  p2 <- Problem(Minimize(y), list(x+y >= 3, y >= 4))
  # g <- partial_optimize(p2, list(y), list(x))
  # p3 <- Problem(Minimize(x+g), list(x >= 5))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
})

test_that("test partial_optimize special constr", {
  x <- Variable(1)
  y <- Variable(1)
  
  # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
  p1 <- Problem(Minimize(x+exp(y)), list(x+y >= 3, y >= 4, x >= 5))
  # result1 <- solve(p1)
  
  # Solve the two-stage problem via partial_optimize
  p2 <- Problem(Minimize(exp(y)), list(x+y >= 3, y >= 4))
  # g <- partial_optimize(p2, list(y), list(x))
  # p3 <- Problem(Minimize(x+g), list(x >= 5))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
})

test_that("test partial_optimize with parameters", {
  x <- Variable(1)
  y <- Variable(1)
  gamma <- Parameter()
  
  # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
  p1 <- Problem(Minimize(x+y), list(x+y >= gamma, y >= 4, x >= 5))
  gamma@value <- 3
  # result1 <- solve(p1)
  
  # Solve the two-stage problem via partial_optimize
  p2 <- Problem(Minimize(y), list(x+y >= gamma, y >= 4))
  # g <- partial_optimize(p2, list(y), list(x))
  # p3 <- Problem(Minimize(x+g), list(x >= 5))
  # result3 <- solve(p3)
  # expect_equal(result1$optimal_value, result3$optimal_value)
})

test_that("test partial_optimize numeric function", {
  x <- Variable(1)
  y <- Variable(1)
  xval <- 4
  
  # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
  p1 <- Problem(Minimize(y), list(xval+y >= 3))
  # result1 <- solve(p1)
  
  # Solve the two-stage problem via partial_optimize
  constr <- list(y >= -100)
  p2 <- Problem(Minimize(y), c(x+y >= 3, constr))
  # g <- partial_optimize(p2, list(y), list(x))
  # x@value <- xval
  # y@value <- 42
  # const[1]@dual_variable@value <- 42
  # result <- g@value
  # expect_equal(result, result1$optimal_value)
  # expect_equal(y@value, 42)
  # expect_equal(constr[1]@dual_value, 42)
  
  # No variables optimized over
  p2 <- Problem(Minimize(y), list(x+y >= 3))
  # g <- partial_optimize(p2, list(), list(x,y))
  # x@value <- xval
  # y@value <- 42
  # p2@constraints[1]@dual_variable@value <- 42
  # result <- g@value
  # expect_equal(result, y@value)
  # expect_equal(y@value, 42)
  # expect_equal(p2@constraints[1]@dual_value, 42)
})

test_that("test partial_optimize stacked", {
  # Minimize the 1-norm in the usual way
  dims <- 3
  x <- Variable(dims)
  t <- Variable(dims)
  p1 <- Problem(Minimize(SumEntries(t)), list(-t <= x, x <= t))
  
  # Minimize the 1-norm via partial_optimize
  # g <- partial_optimize(p1, list(t), list(x))
  # g2 <- partial_optimize(Problem(Minimize(g)), list(x))
  # p2 <- Problem(Minimize(g2))
  # result2 <- solve(p2)
  
  # result1 <- solve(p1)
  # expect_equal(result1$optimal_value, result2$optimal_value)
})

test_that("test the NonNegative Variable class", {
  x <- NonNegative()
  p <- Problem(Minimize(5+x), list(x >= 3))
  # result <- solve(p)
  # expect_equal(result$optimal_value, 8, tolerance = TOL)
  # expect_equal(result$x, 3, tolerance = TOL)
})

test_that("test whether changing an array constant breaks DCP", {
  c <- matrix(c(1, 2), nrow = 2, ncol = 1)
  x@primal_value <- c(1, 1)
  expr <- t(c) %*% Square(x)
  expect_equal(value(expr), 3, tolerance = TOL)
  expect_true(is_dcp(expr))
  
  c[1] <- -1
  expect_equal(value(expr), 3, tolerance = TOL)
  expect_true(is_dcp(expr))
})
