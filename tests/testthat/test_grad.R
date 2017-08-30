TOL <- 1e-6

a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test gradient for AffineProd", {
  value(C) <- rbind(c(1,-2), c(3,4), c(-1,-3))
  value(A) <- rbind(c(3,2), c(-5,1))
  expr <- AffineProd(C, A)
  
  expect_equal(as.matrix(grad(expr)[[as.character(C@id)]]), rbind(c(3,0,0,2,0,0), c(0,3,0,0,2,0), c(0,0,3,0,0,2),
                                                    c(-5,0,0,1,0,0), c(0,-5,0,0,1,0), c(0,0,-5,0,0,1)), tolereance = TOL)
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(1,3,-1,0,0,0), c(-2,4,-3,0,0,0), c(0,0,0,1,3,-1), 
                                                   c(0,0,0,-2,4,-3)), tolerance = TOL)
})

test_that("Test gradient for Pnorm", {
  value(x) <- c(-1,0)
  expr <- Pnorm(x, 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(-1,0)))
  
  value(x) <- c(0,10)
  expr <- Pnorm(x, 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(0,1)))
  
  value(x) <- c(-3,4)
  expr <- Pnorm(x, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(-3.0/5,4.0/5)))
  
  value(x) <- c(-1,2)
  expr <- Pnorm(x, 0.5)
  expect_null(grad(expr)[[as.character(x@id)]])
  
  value(x) <- c(0,0)
  expr <- Pnorm(x, 0.5)
  expect_null(grad(expr)[[as.character(x@id)]])
  
  value(x) <- c(0,0)
  expr <- Pnorm(x, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(0,0)))
  
  value(A) <- rbind(c(2,-2), c(2,2))
  expr <- Pnorm(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0.5,0.5,-0.5,0.5)))
  
  value(A) <- rbind(c(3,-3), c(4,4))
  expr <- Pnorm(A, 2, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0.6,0), c(0.8,0), c(0,-0.6), c(0,0.8)))
  
  value(A) <- rbind(c(3,-4), c(4,3))
  expr <- Pnorm(A, 2, axis = 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0.6,0), c(0,0.8), c(-0.8,0), c(0,0.6)))
  
  value(A) <- rbind(c(3,-4), c(4,3))
  expr <- Pnorm(A, 0.5)
  expect_null(grad(expr)[[as.character(A@id)]])
})

test_that("Test gradient for LogSumExp", {
  value(x) <- c(0,1)
  expr <- LogSumExp(x)
  e <- exp(1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1.0/(1+e), e/(1+e))))
  
  value(A) <- rbind(c(0,1), c(-1,0))
  expr <- LogSumExp(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1.0/(2+e+1.0/e), 1.0/e/(2+e+1.0/e), e/(2+e+1.0/e), 1.0/(2+e+1.0/e))))
  
  value(A) <- rbind(c(0,1), c(-1,0))
  expr <- LogSumExp(A, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), cbind(c(1.0/(1+1.0/e), 1.0/e/(1+1.0/e), 0, 0), c(0, 0, e/(1+e), 1.0/(1+e))))
})

test_that("Test gradient for GeoMean", {
  value(x) <- c(1,2)
  expr <- GeoMean(x)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(sqrt(2)/2, 1.0/2/sqrt(2))))
  
  value(x) <- c(0,2)
  expr <- GeoMean(x)
  expect_null(grad(expr)[[as.character(x@id)]])
  
  value(x) <- c(1,2)
  expr <- GeoMean(x, c(1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,0)))
  
  # No exception for single weight
  value(x) <- c(-1,2)
  expr <- GeoMean(x, c(1,0))
  expect_null(grad(expr)[[as.character(x@id)]])
})

test_that("Test gradient for LambdaMax", {
  value(A) <- rbind(c(2,0), c(0,1))
  expr <- LambdaMax(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,0)))
  
  value(A) <- rbind(c(1,0), c(0,2))
  expr <- LambdaMax(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))
  
  value(A) <- rbind(c(1,0), c(0,1))
  expr <- LambdaMax(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))
})

test_that("Test gradient for MatrixFrac", {
  value(A) <- diag(2)
  value(B) <- diag(2)
  expr <- MatrixFrac(A, B)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(2,0,0,2)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(B@id)]]), matrix(c(-1,0,0,-1)))
  
  value(B) <- matrix(0, nrow = 2, ncol = 2)
  expr <- MatrixFrac(A, B)
  expect_null(grad(expr)[[as.character(A@id)]])
  expect_null(grad(expr)[[as.character(B@id)]])
  
  value(x) <- c(2,3)
  value(A) <- diag(2)
  expr <- MatrixFrac(x, A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(4,6)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(-4,-6,-6,-9)))
})

test_that("Test gradient for NormNuc", {
  value(A) <- rbind(c(10,4), c(4,30))
  expr <- NormNuc(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,1)))
})

test_that("Test gradient for LogDet", {
  value(A) <- 2*diag(rep(1,2))
  expr <- LogDet(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), 1.0/2*diag(2))
  
  mat <- rbind(c(1,2), c(3,5))
  value(A) <- t(mat) %*% mat
  expr <- LogDet(A)
  val <- t(base::solve(value(A)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  
  value(A) <- matrix(0, nrow = 2, ncol = 2)
  expr <- LogDet(A)
  expect_null(grad(expr)[[as.character(A@id)]])
  
  value(A) <- -rbind(c(1,2), c(3,4))
  expr <- LogDet(A)
  expect_null(grad(expr)[[as.character(A@id)]])
})

test_that("Test gradient for QuadOverLin", {
  value(x) <- c(1,2)
  value(a) <- 2
  expr <- QuadOverLin(x, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,2)))
  expect_equal(grad(expr)[[as.character(a@id)]], -1.25)
  
  value(a) <- 0
  expr <- QuadOverLin(x, a)
  expect_null(grad(expr)[[as.character(x@id)]])
  expect_null(grad(expr)[[as.character(a@id)]])
  
  value(A) <- diag(rep(1,2))
  value(a) <- 2
  expr <- QuadOverLin(A, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,1)))
  expect_equal(grad(expr)[[as.character(a@id)]], -0.5)
  
  value(x) <- c(1,2)
  value(a) <- 2
  value(y) <- c(1,2)
  value(a) <- 2
  expr <- QuadOverLin(x, a) + QuadOverLin(y, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,2)))
  expect_equivalent(as.matrix(grad(expr)[[y@id]]), matrix(c(1,2)))
  expect_equal(grad(expr)[[as.character(a@id)]], -2.5)
})

test_that("Test gradient for MaxEntries", {
  value(x) <- c(2,1)
  expr <- MaxEntries(x)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,0)))
  
  value(A) <- rbind(c(1,2), c(4,3))
  expr <- MaxEntries(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,1,0,0)))
  
  value(A) <- rbind(c(1,2), c(4,3))
  expr <- MaxEntries(A, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0,0), c(1,0), c(0,0), c(0,1)))
  
  value(A) <- rbind(c(1,2), c(4,3))
  expr <- MaxEntries(A, axis = 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0,0), c(0,1), c(1,0), c(0,0)))
})

test_that("Test SigmaMax", {
  value(A) <- rbind(c(1,0), c(0,2))
  expr <- SigmaMax(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))
  
  value(A) <- rbind(c(1,0), c(0,1))
  expr <- SigmaMax(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,0)))
})

test_that("Test SumLargest", {
  value(A) <- rbind(c(4,3), c(2,1))
  expr <- SumLargest(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,1,0)))
  
  value(A) <- rbind(c(1,2), c(3,0.5))
  expr <- SumLargest(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,1,1,0)))
})

test_that("Test Abs", {
  value(A) <- rbind(c(1,2), c(-1,0))
  expr <- Abs(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,1,-1,0))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test linearize method", {
  # Affine
  value(x) <- c(1,2)
  expr <- (2*x - 5)[1]
  lin_expr <- linearize(expr)
  value(x) <- c(55,22)
  expr <- (2*x - 5)[1]
  expect_equal(value(lin_expr), value(expr))
  value(x) <- c(-1,-5)
  expr <- (2*x - 5)[1]
  expect_equal(value(lin_expr), value(expr))
  
  # Convex
  expr <- A^2 + 5
  expect_error(linearize(expr))
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- A^2 + 5
  lin_expr <- linearize(expr)
  manual <- value(expr) + 2*Reshape(value(Diag(Vec(A))) * Vec(A - value(A)), 2, 2)
  expect_equal(value(lin_expr), value(expr))
  
  value(A) <- rbind(c(-5,-5), c(8.2,4.4))
  expr <- A^2 + 5
  expect_true(all(value(lin_expr) <= value(expr)))
  expect_equal(value(lin_expr), value(manual))
  
  # Concave
  value(x) <- c(1,2)
  expr <- Log(x)/2
  lin_expr <- linearize(expr)
  manual <- value(expr) + value(Diag(0.5*x^-1))*(x - value(x))
  expect_equal(value(lin_expr), value(expr))
  
  value(x) <- c(3,4.4)
  expr <- Log(x)/2
  expect_true(all(value(lin_expr) >= value(expr)))
  expect_equal(value(lin_expr), value(manual))
})

test_that("Test gradient for Log", {
  value(a) <- 2
  expr <- Log(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/2)
  
  value(a) <- 3
  expr <- Log(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/3)
  
  value(a) <- -1
  expr <- Log(a)
  expect_equal(grad(expr)[[as.character(a@id)]], NA)
  
  value(x) <- c(3,4)
  expr <- Log(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/3,1/4))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(1e-9,4)
  expr <- Log(x)
  expect_equal(grad(expr)[[as.character(x@id)]], NA)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Log(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1, 1/2, 1/3, 1/4))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for Log1p", {
  value(a) <- 2
  expr <- Log1p(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/3)
  
  value(a) <- 3
  expr <- Log1p(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/4)
  
  value(a) <- -1
  expr <- Log1p(a)
  expect_equal(grad(expr)[[as.character(a@id)]], NA)
  
  value(x) <- c(3,4)
  expr <- Log1p(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/4,1/5))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(-1e-9-1,4)
  expr <- Log1p(x)
  expect_equal(grad(expr)[[as.character(x@id)]], NA)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Log1p(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1/2, 1/3, 1/4, 1/5))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for Entr", {
  value(a) <- 2
  expr <- Entr(a)
  expect_equal(grad(expr)[[as.character(a@id)]], -log(2)-1)
  
  value(a) <- 3
  expr <- Entr(a)
  expect_equal(grad(expr)[[as.character(a@id)]], -(log(3)+1))
  
  value(a) <- -1
  expr <- Entr(a)
  expect_null(grad(expr)[[as.character(a@id)]])
  
  value(x) <- c(3,4)
  expr <- Entr(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(-(log(c(3,4)) + 1))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(-1e-9,4)
  expr <- Entr(x)
  expect_null(grad(expr)[[as.character(x@id)]])
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Entr(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(-(log(1:4)+1))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for Exp", {
  value(a) <- 2
  expr <- Exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(2))
  
  value(a) <- 3
  expr <- Exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(3))
  
  value(a) <- -1
  expr <- Exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(-1))
  
  value(x) <- c(3,4)
  expr <- Exp(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4)))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(-1e-9,4)
  expr <- Exp(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4)))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Exp(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for logistic", {
  value(a) <- 2
  expr <- Logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(2)/(1+exp(2)))
  
  value(a) <- 3
  expr <- Logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(3)/(1+exp(3)))
  
  value(a) <- -1
  expr <- Logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(-1)/(1+exp(-1)))
  
  value(x) <- c(3,4)
  expr <- Logistic(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4))/(1+exp(3,4)))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(-1e-9,4)
  expr <- Logistic(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4))/(1+exp(c(-1e-9,4))))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Logistic(A)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4)/(1+exp(1:4)))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for Huber", {
  value(a) <- 2
  expr <- Huber(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 2)
  
  value(a) <- 3
  expr <- Huber(a, M = 2)
  expect_equal(grad(expr)[[as.character(a@id)]], 4)
  
  value(a) <- -1
  expr <- Huber(a, M = 2)
  expect_equal(grad(expr)[[as.character(a@id)]], -2)
  
  value(x) <- c(3,4)
  expr <- Huber(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,2))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(x) <- c(-1e-9,4)
  expr <- Huber(x)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,2))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- Huber(A, M = 3)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,4,6,6))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for KLDiv", {
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- KLDiv(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], log(2/4))
  expect_equal(grad(expr)[[as.character(b@id)]], 1-(2/4))
  
  value(a) <- 3
  value(b) <- 0
  expr <- KLDiv(a, b)
  expect_null(grad(expr)[[as.character(a@id)]])
  expect_null(grad(expr)[[as.character(b@id)]])
  
  value(a) <- -1
  value(b) <- 2
  expr <- KLDiv(a, b)
  expect_null(grad(expr)[[as.character(a@id)]])
  expect_null(grad(expr)[[as.character(b@id)]])
  
  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,8)
  expr <- KLDiv(x, y)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(log(c(3,4)) - log(c(5,8)))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <-matrix(0, nrow = 2, ncol = 2) + diag(c(1-3/5,1-4/8))
  expect_equal(as.matrix(grad(expr)[[as.character(y@id)]]), val)
  
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,2)
  expr <- KLDiv(x, y)
  expect_equal(grad(expr)[[as.character(x@id)]], NA)
  expect_equal(grad(expr)[[y@id]], NA)
  
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3.5,2.3))
  expr <- KLDiv(A, B)
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(log(div))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(1-div)
  expect_equal(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for MaxElemwise", {
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- MaxElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)
  
  value(a) <- 3
  value(b) <- 0
  expr <- MaxElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)
  
  value(a) <- -1
  value(b) <- 2
  expr <- MaxElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)
  
  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  expr <- MaxElemwise(x, y)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  expr <- MaxElemwise(x, y)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3,2.3))
  expr <- MaxElemwise(A, B)
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,1,1))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,0,0))
  expect_equal(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for MinElemwise", {
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- MinElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)
  
  value(a) <- 3
  value(b) <- 0
  expr <- MinElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)
  
  value(a) <- -1
  value(b) <- 2
  expr <- MinElemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)
  
  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  expr <- MinElemwise(x, y)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  expr <- MinElemwise(x, y)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3,2.3))
  expr <- MinElemwise(A, B)
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,1,0))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,0,1))
  expect_equal(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for Power", {
  value(a) <- 2
  expr <- Sqrt(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5/sqrt(2))
  
  value(a) <- 3
  expr <- Sqrt(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5/sqrt(3))
  
  value(a) <- -1
  expr <- Sqrt(a)
  expect_null(grad(expr)[[as.character(a@id)]])
  
  value(x) <- c(3,4)
  expr <- x^3
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), rbind(c(27,0), c(0,48)))
  
  value(x) <- c(-1e-9,4)
  expr <- x^3
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), rbind(c(0,0), c(0,48)))
  
  value(A) <- rbind(c(1,-2), c(3,4))
  expr <- A^2
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(2,-4,6,8))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  
  # Constant
  expr <- a^0
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  
  expr <- x^0
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(0, nrow = 2, ncol = 2))
})

test_that("Test grad for partial minimization/maximization problems", {
  # for(obj in list(Minimize(a^-1), Maximize(Entr(a)))) {
  #   prob <- Problem(obj, list(x + a >= c(5,8)))
  #   
  #   # Optimize over nothing
  #   expr <- partial_optimize(prob, dont_opt_vars = list(x, a))
  #   value(a) <- NA
  #   value(x) <- NA
  #   grad <- grad(expr)
  #   expect_null(grad[[as.character(a@id)]])
  #   expect_null(grad[[as.character(x@id)]])
  #   
  #   # Outside domain
  #   value(a) <- 1.0
  #   value(x) <- c(5,5)
  #   grad <- grad(expr)
  #   expect_null(grad[[as.character(a@id)]])
  #   expect_null(grad[[as.character(x@id)]])
  #   
  #   value(a) <- 1
  #   value(x) <- c(10,10)
  #   grad <- grad(expr)
  #   expect_equal(grad[[as.character(a@id)]], grad(obj@.args[[1]])[[as.character(a@id)]])
  #   expect_equal(as.matrix(grad[[as.character(x@id)]]), rep(0,4))
  #   
  #   # Optimize over x
  #   expr <- partial_optimize(prob, opt_vars = list(x))
  #   value(a) <- 1
  #   grad <- grad(expr)
  #   expect_equal(grad[[as.character(a@id)]], grad(obj@.args[[1]])[[as.character(a@id)]] + 0)
  #   
  #   # Optimize over a
  #   fix_prob <- Problem(obj, list(x + a >= c(5,8), x == 0))
  #   result <- solve(fix_prob)
  #   dual_val <- fix_prob@constraints[1]@dual_variable@value
  #   expr <- partial_optimize(prob, opt_vars = list(a))
  #   value(x) <- c(0,0)
  #   grad <- grad(expr)
  #   expect_equal(as.matrix(grad[[as.character(x@id)]]), dual_val)
  #   
  #   # Optimize over x and a
  #   expr <- partial_optimize(prob, opt_vars = list(x, a))
  #   grad <- grad(expr)
  #   expect_equal(grad, list())
  # }
})

test_that("Test grad for affine atoms", {
  value(a) <- 2
  expr <- -a
  expect_equal(grad(expr)[[as.character(a@id)]], -1)
  
  value(a) <- 2
  expr <- 2*a
  expect_equal(grad(expr)[[as.character(a@id)]], 2)
  
  value(a) <- 2
  expr <- a/2
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5)
  
  value(x) <- c(3,4)
  expr <- -x
  val <- matrix(0, nrow = 2, ncol = 2) - diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- -A
  val <- matrix(0, nrow = 4, ncol = 4) - diag(rep(1,4))
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  
  value(A) <- rbind(c(1,2), c(3,4))
  expr <- A[1,2]
  val <- matrix(0, nrow = 4, ncol = 1)
  val[3] <- 1
  expect_equal(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  
  z <- Variable(3)
  value(x) <- c(1,2)
  value(z) <- c(1,2,3)
  expr <- VStack(x, z)
  val <- matrix(0, nrow = 2, ncol = 5)
  val[,1:3] <- diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  
  val <- matrix(0, nrow = 3, ncol = 5)
  val[,3:ncol(val)] <- diag(rep(1,3))
  expect_equal(as.matrix(grad(expr)[[as.character(z@id)]]), val)
})
