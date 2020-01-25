context("test-benchmarks")

if(!require(microbenchmark))
  install.packages("microbenchmark")

perform <- CVXR:::perform
ConeMatrixStuffing <- CVXR:::ConeMatrixStuffing

test_that("test diffcp SDP example", {
  randn_symm <- function(n) {
    A <- matrix(rnorm(n*n), nrow = n, ncol = n)
    return((A + t(A))/2)
  }

  randn_psd <- function(n) {
    A <- 1.0 / (10 * matrix(rnorm(n*n), nrow = n, ncol = n))
    return(A %*% t(A))
  }

  n <- 100
  p <- 100
  C <- randn_psd(n)
  As <- lapply(1:p, function(i) { randn_symm(n) })
  Bs <- matrix(rnorm(p), nrow = p)

  diffcp_sdp <- function() {
    X <- Variable(n, n, PSD = TRUE)
    objective <- matrix_trace(C %*% X)
    constraints <- lapply(1:p, function(i) { As[[i]] %*% X == Bs[i] })
    problem <- Problem(Minimize(objective), constraints)
    get_problem_data(problem, "SCS")
  }
  microbenchmark(diffcp_sdp(), times = 1L)
})

test_that("test TV inpainting", {
    ## 512 is too big!
    ## rows <- 512
    ## cols <- 512
    ##
    rows <- 64
    cols <- 64
    colors <- 3

  Uorig <- array(rnorm(rows*cols*colors), dim = c(rows, cols, colors))
  known <- array(0, dim = c(rows, cols, colors))
  for(i in seq_len(rows)) {
    for(j in seq_len(cols)) {
      if(runif(1) > 0.7) {
        for(k in seq_len(colors))
          known[i,j,k] <- 1
      }
    }
  }

  tv_inpainting <- function() {
    Ucorr <- known*Uorig
    variables <- list()
    constraints <- list()

    for(i in seq_len(colors)) {
      U <- Variable(rows, cols)
      variables <- c(variables, U)
      constraints <- c(constraints, known[,,i] * U == known[,,i] * Ucorr[,,i])
    }
    obj <- do.call(tv, variables)
    prob <- Problem(Minimize(obj), constraints)
    get_problem_data(prob, "SCS")
  }
  microbenchmark(tv_inpainting(), times = 3L)
})

test_that("test least squares", {
  m <- 20
  n <- 15
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m), nrow = m)

  least_squares <- function() {
    x <- Variable(n)
    cost <- sum_squares(A %*% x - b)
    get_problem_data(Problem(Minimize(cost)), "OSQP")
  }
  microbenchmark(least_squares(), times = 1L)
})

test_that("test qp", {
  m <- 15
  n <- 10
  p <- 5
  P <- matrix(rnorm(n*n), nrow = n, ncol = n)
  P <- t(P) %*% P
  q <- matrix(rnorm(n), nrow = n)
  G <- matrix(rnorm(m*n), nrow = m, ncol = n)
  h <- G %*% matrix(rnorm(n), nrow = n)
  A <- matrix(rnorm(p*n), nrow = p, ncol = n)
  b <- matrix(rnorm(p), nrow = p)

  qp <- function() {
    x <- Variable(n)
    get_problem_data(Problem(Minimize((1/2)*quad_form(x, P) + t(q) %*% x),
                             list(G %*% x <= h, A %*% x == b)), "OSQP")
  }
  microbenchmark(qp(), times = 1L)
})

test_that("test stuffing perf with many constraints", {
  m <- 2000
  n <- 2000
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  C <- matrix(runif(floor(m/2)), nrow = floor(m/2))
  b <- matrix(rnorm(m), nrow = m)

  x <- Variable(n)
  cost <- sum(A %*% x)

  constraints <- lapply(1:floor(m/2), function(i) { C[i] * x[i] <= b[i] })
  constraints <- c(constraints, lapply(1:floor(m/2), function(i) { C[i] * x[floor(m/2) + i] == b[floor(m/2) + i] }))

  p <- Problem(Minimize(cost), constraints)

  stuff <- function(mat) {
    perform(ConeMatrixStuffing(), mat)
  }

  microbenchmark(stuff(p), times = 1L)
})
