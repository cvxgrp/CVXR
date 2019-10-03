# TODO_NARAS_12: Should we retire this test since it no longer exists in CVXPY?
context("test-g03-ls")
TOL <- 1e-6

test_that("Test regression", {
    ## Set the random seed to get consistent data
    skip_on_cran()
    set.seed(1)

    ## Number of examples to use
    n <- 100

    ## Specify the true value of the variable
    true_coeffs <- matrix(c(2, -2, 0.5), nrow = 3, ncol = 1)

    ## Generate data
    x_data <- matrix(stats::runif(n) * 5, nrow = n, ncol = 1)
    x_data_expanded <- cbind(x_data, x_data^2, x_data^3)
    y_data <- x_data_expanded %*% true_coeffs + 0.5 * matrix(stats::runif(n, 1), nrow = n, ncol = 1)

    slope <- Variable()
    offset <- Variable()
    line <- offset + x_data * slope
    residuals <- line - y_data
    fit_error <- sum_squares(residuals)
    result <- solve(Problem(Minimize(fit_error), list()))
    ## result <- solve(Problem(Minimize(fit_error), list()), solver = "LS")
    ## optval <- result$value
    ## expect_equal(optval, 1171.60037715, tolerance = TOL)

    quadratic_coeff <- Variable()
    slope <- Variable()
    offset <- Variable()
    quadratic <- offset + x_data * slope + quadratic_coeff * x_data^2
    residuals <- quadratic - y_data
    fit_error <- sum_squares(residuals)
    ## result <- solve(Problem(Minimize(fit_error), list()), solver = "LS")
    result2 <- solve(Problem(Minimize(fit_error), list()), solver = "ECOS")
    ## optval <- result$value
    optval2 <- result2$value
    ## expect_equal(optval, 139.225650756, tolerance = TOL)
    expect_equal(optval2, 109.639449308, tolerance = TOL)
})

test_that("Test control", {
    ## Some constraints on our motion
    ## The object should start from the origin and end at rest
    skip_on_cran()
    initial_velocity <- matrix(c(-20, 100), nrow = 2, ncol = 1)
    final_position <- matrix(c(100, 100), nrow = 2, ncol = 1)

    Tnum <- 5  ## The number of timesteps
    h <- 0.1   ## The time between time intervals
    mass <- 1  ## Mass of object
    drag <- 0.1  ## Drag on object
    g <- matrix(c(0, -9.8), nrow = 2, ncol = 1)   ## Gravity on object

    ## Declare the variables we need
    position <- Variable(2, Tnum)
    velocity <- Variable(2, Tnum)
    force <- Variable(2, Tnum-1)

    ## Create a problem instance
    mu <- 1
    constraints <- list()

    ## Add constraints on our variables
    for(i in 1:(Tnum-1)) {
        constraints <- c(constraints, position[,i+1] == position[,i] + h * velocity[,i])
        acceleration <- force[,i]/mass + g - drag * velocity[,i]
        constraints <- c(constraints, velocity[,i+1] == velocity[,i] + h * acceleration)
    }

    ## Add position constraints
    constraints <- c(constraints, position[,1] == 0)
    constraints <- c(constraints, position[,Tnum] == final_position)

    ## Add velocity constraints
    constraints <- c(constraints, velocity[,1] == initial_velocity)
    constraints <- c(constraints, velocity[,Tnum] == 0)

    ## Solve the problem
    result <- solve(Problem(Minimize(sum_squares(force)), constraints))   ## TODO: Remove once LS solver is implemented.
    ## result <- solve(Problem(Minimize(sum_squares(force)), constraints), solver = "LS")
    ## optval <- result$value
    ## expect_equal(optval, 17859.0, tolerance = 1)
})

test_that("Test sparse system", {
    skip_on_cran()
    m <- 1000
    n <- 800
    r <- 700

    set.seed(1)
    density <- 0.2
    A <- Matrix::rsparsematrix(m, n, density, rand.x = stats::runif)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)
    G <- Matrix::rsparsematrix(r, n, density, rand.x = stats::runif)
    h <- matrix(stats::rnorm(r), nrow = r, ncol = 1)

    x <- Variable(n)
    result <- solve(Problem(Minimize(sum_squares(A %*% x - b)), list(G %*% x == h)))
    ## result <- solve(Problem(Minimize(sum_squares(A %*% x - b)), list(G %*% x == h)), solver = "LS")
    optval <- result$value
    expect_equal(optval, 6542.100424497, tolerance = TOL)
})

test_that("Test equivalent forms", {
    skip_on_cran()
    m <- 100
    n <- 80
    r <- 70

    set.seed(1)
    A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)
    G <- matrix(stats::rnorm(r*n), nrow = r, ncol = n)
    h <- matrix(stats::rnorm(r), nrow = r, ncol = 1)

    ## ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2
    P <- t(A) %*% A
    q <- -2 * t(A) %*% b
    r <- t(b) %*% b

    Pinv <- base::solve(P)

    x <- Variable(n)

    obj1 <- sum_squares(A %*% x - b)
    obj2 <- sum_entries((A %*% x - b)^2)
    obj3 <- quad_form(x, P) + t(q) %*% x + r
    obj4 <- matrix_frac(x, Pinv) + t(q) %*% x + r

    cons <- list(G %*% x == h)

    v1 <- solve(Problem(Minimize(obj1), cons), solver = "ECOS")$value
    v2 <- solve(Problem(Minimize(obj2), cons), solver = "ECOS")$value
    v3 <- solve(Problem(Minimize(obj3), cons), solver = "ECOS")$value
    ## v1 <- solve(Problem(Minimize(obj1), cons), solver = "LS")$value
    ## v2 <- solve(Problem(Minimize(obj2), cons), solver = "LS")$value
    ## v3 <- solve(Problem(Minimize(obj3), cons), solver = "LS")$value
    ## v4 <- solve(Problem(Minimize(obj4), cons), solver = "LS")$value

    expect_equal(v1, 622.204576157, tolerance = TOL)
    expect_equal(v2, 622.204576157, tolerance = TOL)
    expect_equal(v3, 622.204576157, tolerance = TOL)
    ## expect_equal(v4, 622.204576157, tolerance = TOL)
})

test_that("Test smooth ridge", {
    skip_on_cran()
    set.seed(1)
    n <- 500
    k <- 50
    delta <- 1
    eta <- 1

    A <- matrix(stats::runif(k*n), nrow = k, ncol = n)
    b <- matrix(stats::runif(k), nrow = k, ncol = 1)
    x <- Variable(n)
    obj <- sum_squares(A %*% x - b) + delta*sum_squares(x[1:(n-1)]-x[2:n]) + eta*sum_squares(x)
    optval <- solve(Problem(Minimize(obj)))$value
    ## optval <- solve(Problem(Minimize(obj), list()), solver = "LS")$value
    expect_equal(optval, 0.208432916, tolerance = TOL)
})
