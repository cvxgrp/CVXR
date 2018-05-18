## These tests were moved here because scs causes errors on R-devel on travis
## Also because they go over the time limit

TOL <- 1e-6


test_that("Test examples from the README", {
    ## Problem data
    skip_on_cran()
    m <- 30
    n <- 20
    A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)

    ## Construct the problem
    x <- Variable(n)
    objective <- Minimize(sum((A %*% x - b)^2))
    constraints <- list(x >= 0, x <= 1)
    p <- Problem(objective, constraints)

    ## The optimal objective is returned by solve(p)
    result <- solve(p)
    ## The optimal value for x is stored in result$getValue(x)
    print(result$getValue(x))
    ## The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraints[[1]])
    print(result$getDualValue(constraints[[1]]))

###########################################
    ## Scalar variable
    a <- Variable()

    ## Column vector variable of length 5
    x <- Variable(5)

    ## Matrix variable with 4 rows and 7 columns
    A <- Variable(4, 7)

###########################################
    ## Positive scalar parameter
    m <- Parameter(sign = "positive")

    ## Column vector parameter with unknown sign (by default)
    c <- Parameter(5)

    ## Matrix parameter with negative entries
    G <- Parameter(4, 7, sign = "negative")

    ## Assigns a constant value to G
    value(G) <- -matrix(1, nrow = 4, ncol = 7)

    expect_error(value(G) <- matrix(1, nrow = 4, ncol = 7))

###########################################
    a <- Variable()
    x <- Variable(5)

    ## expr is an Expression object after each assignment
    expr <- 2*x
    expr <- expr - a
    expr <- sum(expr) + norm2(x)

###########################################
    ## Problem data
    n <- 10
    m <- 5
    A <- matrix(stats::rnorm(n*m), nrow = n, ncol = m)
    b <- matrix(stats::rnorm(n), nrow = n, ncol = 1)
    gamma <- Parameter(sign = "positive")

    ## Construct the problem
    x <- Variable(m)
    loss <- sum((A %*% x - b)^2)

    get_x <- function(gamma) {
        objective <- Minimize(loss + gamma*norm1(x))
        p <- Problem(objective)

        result <- solve(p)
        result$getValue(x)
    }

    gammas <- 10^seq(-1, 2, length.out = 2)

    ## Serial computation
    x_values <- sapply(gammas, get_x)

###########################################
    n <- 10

    mu <- matrix(stats::rnorm(n), nrow = 1, ncol = n)
    sigma <- matrix(stats::rnorm(n^2), nrow = n, ncol = n)
    sigma <- t(sigma) %*% sigma
    gamma <- Parameter(sign = "positive")
    value(gamma) <- 1
    x <- Variable(n)

    ## Constants:
    ## mu is the vector of expected returns.
    ## sigma is the covariance matrix.
    ## gamma is a Parameter that trades off risk and return.

    ## Variables:
    ## x i s a vector of stock holdings as fractions of total assets.

    expected_return <- mu %*% x
    risk <- quad_form(x, sigma)

    gamma <- 1
    objective <- Maximize(expected_return - gamma*risk)
    p <- Problem(objective, list(sum(x) == 1))
    result <- solve(p)

    ## The optimal expected return
    print(result$getValue(expected_return))

    ## The optimal risk
    print(result$getValue(risk))

###########################################
    N <- 50
    M <- 40
    n <- 10
    data1 <- lapply(1:N, function(i) { list(1, matrix(stats::rnorm(n, mean = 1.0, sd = 2.0))) })
    data2 <- lapply(1:M, function(i) { list(-1, matrix(stats::rnorm(n, mean = -1.0, sd = 2.0))) })
    data <- c(data1, data2)

    ## Construct problem
    gamma <- Parameter(sign = "positive")
    value(gamma) <- 0.1
    a <- Variable(n)
    b <- Variable()

    slack <- lapply(data, function(x) {
        label <- x[[1]]
        sample <- x[[2]]
        pos(1 - label*(t(sample) %*% a - b))
    })
    gamma <- 0.1
    objective <- Minimize(norm2(a) + gamma*Reduce("+", slack))
    p <- Problem(objective)
    result <- solve(p)

    ## Count misclassifications
    errors <- 0
    for(v in data) {
        label <- v[[1]]
        sample <- v[[2]]
        if(label * result$getValue(t(sample) %*% a - b) < 0)
            errors <- errors + 1
    }

    print(errors)
    print(result$getValue(a))
    print(result$getValue(b))
})

test_that("Test image in-painting", {
    skip_on_cran()
    set.seed(1)
    rows <- 100
    cols <- 100

    ## Load the images and convert to arrays
    Uorig <- matrix(sample(0:255, size = rows * cols, replace = TRUE), nrow = rows, ncol = cols)

    ## Known is 1 if the pixel is known, 0 if the pixel was corrupted
    ## Known <- matrix(0, nrow = rows, ncol = cols)

    ## for(i in 1:rows) {
    ##     for(j in 1:cols) {
    ##         if(stats::runif(1) > 0.7)
    ##             Known[i,j] <- 1
    ##     }
    ## }

    Known <- sapply(seq_len(cols), function(x) as.numeric(stats::runif(n = rows) > 0.7))

    Ucorr <- Known %*% Uorig

    ## Recover the original image using total variation in-painting
    U <- Variable(rows, cols)
    obj <- Minimize(tv(U))
    constraints <- list(Known * U == Known * Ucorr)
    prob <- Problem(obj, constraints)
    res <- solve(prob, solver = "SCS") ## Does testing on travis cause SCS to bomb?
    ## It certainly gives warnings, but that is a known problem with SCS
})

test_that("Test the log_sum_exp function", {
    skip_on_cran()
    set.seed(1)
    m <- 5
    n <- 2
    X <- matrix(1, nrow = m, ncol = n)
    w <- Variable(n)

    expr2 <- lapply(1:m, function(i) { log_sum_exp(vstack(0, X[i,] %*% w)) })
    expr3 <- Reduce("+", expr2)
    obj <- Minimize(expr3)
    p <- Problem(obj)
    res <- solve(p, solver = "SCS", max_iters = 1)
})
